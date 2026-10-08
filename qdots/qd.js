// Hot-injection synthesis of CdSe quantum dots.
//
// The dot population is a set of size bins {r (nm, radius), n (number)}.
// Monomer concentration C is in units of the amount added by one injection.
//
// Growth (diffusion limited, with Gibbs-Thomson solubility):
//   dr/dt = K(T) * (C - Cs(T) * exp(alpha / r)) / r
// Dots smaller than the critical radius alpha / ln(C/Cs) dissolve and feed
// the bigger ones: Ostwald ripening. At high supersaturation the 1/r factor
// lets small dots catch up: size focusing.
//
// Nucleation (classical): J = A * exp(-B / ln^2 S), S = C / Cs.
// Emission: Yu et al. (2003) empirical sizing curve plus a Stokes shift.
(function (root) {
  const ALPHA = 1.2;            // capillary length (nm)
  const GAMMA = 0.012;           // monomer used per nm^3 of dot volume (per unit n)
  const T_REF = 573;             // 300 C

  const Cs = T => 0.03 * Math.exp((T - T_REF) / 38);
  const K  = T => 0.05 * Math.exp((T - T_REF) / 45);

  // diameter (nm) for a first absorption peak at lambda (nm)
  const yuD = l => 1.6122e-9 * l ** 4 - 2.6575e-6 * l ** 3 + 1.6242e-3 * l ** 2 - 0.4277 * l + 41.57;
  // inverse by bisection, valid for D ~ 1.7 to 9 nm
  function absPeak(D) {
    let lo = 420, hi = 670;
    if (D <= yuD(lo)) return lo; if (D >= yuD(hi)) return hi;
    for (let i = 0; i < 40; i++) { const m = (lo + hi) / 2; if (yuD(m) < D) lo = m; else hi = m; }
    return (lo + hi) / 2;
  }
  const STOKES = 12;
  const emission = r => absPeak(2 * r) + STOKES;

  function makeFlask() {
    return { T: 573, setT: 573, C: 0, bins: [], t: 0, injections: 0, burst: 0 };
  }

  // Injection: the mixing of cold precursor gives a short burst during which
// nucleation is easy ("burst nucleation", LaMer). Outside the burst only a very
// high supersaturation nucleates new dots.
function inject(f, amount = 1) { f.C += amount; f.injections++; f.burst = BURST; f.T -= 12 * amount; }
  // Slow dropwise addition: monomer arrives gradually, no mixing burst.
  const DRIP_RATE = 0.02;   // per second
  const BURST = 1.5;
  const SPREAD = 0.12;     // dot-to-dot variation in surface reactivity
  function gauss() { let u = 0; while (!u) u = Math.random(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * Math.random()); }

  function totalN(f) { let s = 0; for (const b of f.bins) s += b.n; return s; }

  function step(f, dt) {
    // heater: temperature relaxes towards the set point (~10 s time constant);
    // an injection of cold precursor knocks it down a little (handled by caller)
    f.T += (f.setT - f.T) * Math.min(1, dt / 10);
    if (f.drip) { f.C += DRIP_RATE * dt; f.dripped = (f.dripped || 0) + DRIP_RATE * dt; }
    const cs = Cs(f.T), k = K(f.T);
    // substep for stability
    const nsub = Math.max(1, Math.ceil(dt / 0.05));
    const h = dt / nsub;
    for (let s = 0; s < nsub; s++) {
      const S = f.C / cs;
      // nucleation
      if (S > 1.5) {
        const l = Math.log(S);
        const J = (f.burst > 0 ? 1 : 0.003) * 5 * Math.exp(-26 / (l * l));
        if (J * h > 1e-3) {
          const rc = ALPHA / l;
          f.bins.push({ r: rc * 1.15, n: J * h, k: Math.max(0.5, 1 + SPREAD * gauss()) });
          f.C -= GAMMA * J * h * (4 / 3) * Math.PI * (rc * 1.15) ** 3;
        }
      }
      // growth / dissolution
      let used = 0;
      for (const b of f.bins) {
        const dr = k * b.k * (f.C - cs * Math.exp(ALPHA / b.r)) / b.r * h;
        const r1 = Math.max(0.05, b.r + dr);
        used += b.n * (4 / 3) * Math.PI * (r1 ** 3 - b.r ** 3);
        b.r = r1;
      }
      f.C = Math.max(0, f.C - GAMMA * used);
      // dissolved dots vanish, returning nothing more (already counted)
      for (let i = f.bins.length - 1; i >= 0; i--) if (f.bins[i].r <= 0.06) f.bins.splice(i, 1);
      f.t += h; if (f.burst > 0) f.burst -= h;
    }
    mergeBins(f, 160);
  }

  function mergeBins(f, max) {
    if (f.bins.length <= max) return;
    f.bins.sort((a, b) => a.r - b.r);
    while (f.bins.length > max) {
      let best = 0, gap = Infinity;
      for (let i = 0; i < f.bins.length - 1; i++) {
        const g = f.bins[i + 1].r - f.bins[i].r;
        if (g < gap) { gap = g; best = i; }
      }
      const a = f.bins[best], b = f.bins[best + 1], n = a.n + b.n;
      // conserve volume
      const r = Math.cbrt((a.n * a.r ** 3 + b.n * b.r ** 3) / n);
      f.bins.splice(best, 2, { r, n, k: (a.n * a.k + b.n * b.k) / n });
    }
  }

  // Ensemble emission spectrum on a grid; each dot has a homogeneous Gaussian line.
  const HOMOG_FWHM = 22;
  function spectrum(f, l0 = 420, l1 = 700, dl = 1) {
    const m = Math.round((l1 - l0) / dl) + 1;
    const y = new Float64Array(m);
    const sig = HOMOG_FWHM / 2.3548;
    for (const b of f.bins) {
      if (b.r < 0.7) continue;                    // tiny clusters barely emit
      const lc = emission(b.r);
      const wt = b.n * Math.min(1, (b.r - 0.7) / 0.4);
      for (let i = 0; i < m; i++) { const u = (l0 + i * dl - lc) / sig; if (Math.abs(u) < 4) y[i] += wt * Math.exp(-0.5 * u * u); }
    }
    let peak = 0, ip = 0; for (let i = 0; i < m; i++) if (y[i] > peak) { peak = y[i]; ip = i; }
    let fwhm = 0;
    if (peak > 0) {
      let a = ip, b = ip;
      while (a > 0 && y[a] > peak / 2) a--;
      while (b < m - 1 && y[b] > peak / 2) b++;
      fwhm = (b - a) * dl;
    }
    return { l0, dl, y, peak, lambda: l0 + ip * dl, fwhm };
  }

  function sizeStats(f) {
    let n = 0, s = 0, s2 = 0;
    for (const b of f.bins) { if (b.r < 0.7) continue; const d = 2 * b.r; n += b.n; s += b.n * d; s2 += b.n * d * d; }
    if (!n) return { mean: 0, rsd: 0, n: 0 };
    const mean = s / n;
    return { mean, rsd: Math.sqrt(Math.max(0, s2 / n - mean * mean)) / mean, n };
  }

  root.QD = { makeFlask, inject, step, spectrum, sizeStats, emission, Cs, K, totalN };
  if (typeof module !== 'undefined') module.exports = root.QD;
})(typeof window !== 'undefined' ? window : globalThis);
