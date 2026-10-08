// Laser micromachining of a metal surface, in cross-section.
// The surface is a height profile h(x) in micrometres (0 = original surface,
// negative = material removed, positive = burr / recast).
//
// Each pulse has a Gaussian fluence profile F(r) = F0 exp(-2 r^2 / w0^2).
// Ablation depth per pulse follows the logarithmic law
//   dz = delta * ln(F_eff / F_th)   for F_eff > F_th, else 0,
// where F_eff = F * cos(theta) accounts for the beam hitting a sloped wall.
// Fluences are in units of the threshold F_th.
(function (root) {
  const LASERS = {
    fs: { name: 'Femtosecond', wavelength: '1030 nm', Fth: 0.2, delta: 0.04, burr: 0, haz: 0 },
    ns: { name: 'Nanosecond', wavelength: '1064 nm', Fth: 2.0, delta: 0.30, burr: 0.45, haz: 1 },
  };

  function makeWork(widthUm, dx, depthUm) {
    const n = Math.round(widthUm / dx) + 1;
    return { dx, n, width: widthUm, floor: -depthUm, h: new Float64Array(n), heat: new Float64Array(n), pulses: 0 };
  }

  function reset(w) { w.h.fill(0); w.heat.fill(0); w.pulses = 0; }

  // Fire one pulse centred at xc (um), spot radius w0 (um), peak fluence F0 (x threshold).
  // Returns the volume removed (um^2 per unit length) so callers can size the plume.
  function pulse(w, laser, xc, w0, F0) {
    const { h, dx, n } = w;
    const L = LASERS[laser];
    const dz = new Float64Array(n);
    let removed = 0;
    const i0 = Math.max(1, Math.floor((xc - 2 * w0) / dx));
    const i1 = Math.min(n - 2, Math.ceil((xc + 2 * w0) / dx));
    for (let i = i0; i <= i1; i++) {
      const r = i * dx - xc;
      const F = F0 * Math.exp(-2 * r * r / (w0 * w0));
      if (F <= 1) continue;
      // local wall slope from a ~1 um stencil (a 3-point stencil lets
      // neighbouring points decouple and grow spikes)
      const a = Math.max(0, i - 5), b = Math.min(n - 1, i + 5);
      const s = (h[b] - h[a]) / ((b - a) * dx);
      const Feff = F / Math.sqrt(1 + s * s);
      if (Feff > 1) dz[i] = L.delta * Math.log(Feff);
    }
    // smooth the removal over ~0.5 um: the real surface can't hold features
    // much finer than the beam's optical resolution
    const sm = new Float64Array(n);
    for (let i = i0; i <= i1; i++) {
      let acc = 0, wsum = 0;
      for (let j = -3; j <= 3; j++) {
        const q = i + j; if (q < 0 || q >= n) continue;
        const wt = 4 - Math.abs(j); acc += wt * dz[q]; wsum += wt;
      }
      sm[i] = Math.min(acc / wsum, h[i] - w.floor);
      removed += sm[i] * dx;
    }
    for (let i = i0; i <= i1; i++) h[i] -= sm[i];
    // surface relaxation where material is being removed (melt flow and
    // redeposition smooth sub-micron roughness); untouched walls are left alone
    for (let it = 0; it < 2; it++) {
      const hc = h.slice(i0 - 1, i1 + 2);
      for (let i = i0; i <= i1; i++) {
        const c = 0.24 * Math.min(1, sm[i] / 0.03);
        if (c > 0) h[i] = hc[i - i0 + 1] + c * (hc[i - i0] - 2 * hc[i - i0 + 1] + hc[i - i0 + 2]);
      }
    }

    if (L.burr > 0 && removed > 0 && F0 > 1) {
      // melt is pushed out and resolidifies as a rim just outside the crater edge;
      // a hotter pulse throws a larger fraction of it out of the hole
      const edge = w0 * Math.sqrt(Math.log(F0) / 2);
      const sig = 0.35 * w0;
      const frac = L.burr * Math.min(1, 0.4 + 0.08 * F0);
      const share = (frac * removed) / 2;
      const norm = share / (sig * Math.sqrt(2 * Math.PI));
      for (const side of [-1, 1]) {
        const xr = xc + side * (edge + 0.5 * sig);
        const j0 = Math.max(0, Math.floor((xr - 3 * sig) / dx));
        const j1 = Math.min(n - 1, Math.ceil((xr + 3 * sig) / dx));
        for (let j = j0; j <= j1; j++) {
          const u = (j * dx - xr) / sig;
          h[j] += norm * Math.exp(-0.5 * u * u) * dx;
        }
      }
    }
    if (L.haz > 0) {
      const reach = 3 * w0;
      for (let i = Math.max(0, Math.floor((xc - reach) / dx)); i <= Math.min(n - 1, Math.ceil((xc + reach) / dx)); i++) {
        const r = i * dx - xc;
        w.heat[i] = Math.min(1, w.heat[i] + 0.004 * F0 * Math.exp(-r * r / (reach * reach / 2)));
      }
    }
    w.pulses++;
    return removed;
  }

  // Crater diameter for a single pulse on a flat surface (the classic D^2 = 2 w0^2 ln(F0/Fth)).
  function craterDiameter(w0, F0) { return F0 > 1 ? 2 * w0 * Math.sqrt(Math.log(F0) / 2) : 0; }

  // Tolerance band around the target outline (um): real laser-cut walls taper.
  function tolerance(tw, td) { return Math.max(0.8, Math.min(2, 0.12 * Math.min(tw, 2 * td))); }

  // Score a profile against a rectangular target channel centred at xc.
  // Anything inside the tolerance band is fine; error is the area outside it.
  function score(w, xc, tw, td) {
    const { h, dx, n } = w;
    const tol = tolerance(tw, td);
    let under = 0, over = 0, burr = 0;
    for (let i = 0; i < n; i++) {
      const ax = Math.abs(i * dx - xc);
      const got = Math.max(0, -h[i]);
      let lo, hi;
      if (ax <= tw / 2 - tol) { lo = td - tol; hi = td + tol; }   // channel floor
      else if (ax <= tw / 2 + tol) { lo = 0; hi = td + tol; }     // wall band
      else { lo = 0; hi = 0.25; }                                  // untouched surface
      if (got < lo) under += (lo - got) * dx;
      if (got > hi) over += (got - hi) * dx;
      if (h[i] > 0.1) burr += h[i] * dx;
    }
    const area = tw * td;
    const err = under + over + 3 * burr;
    const pct = Math.max(0, Math.round(100 * (1 - err / area)));
    return { pct, under, over, burr, tol };
  }

  root.Laser = { LASERS, makeWork, reset, pulse, craterDiameter, score, tolerance };
  if (typeof module !== 'undefined') module.exports = root.Laser;
})(typeof window !== 'undefined' ? window : globalThis);
