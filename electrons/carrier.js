// Free Path: electron and hole transport in an n-type silicon bar.
// Units: cm, s, V, K unless noted. Works in the browser (window.Carrier) and in Node.
(function (root) {
  const q = 1.602e-19, kB = 8.617e-5;          // C, eV/K
  const kJ = 1.381e-23, m0 = 9.109e-31;        // J/K, kg
  const ME = 0.26, MH = 0.36;                  // conductivity effective masses (× m0)

  // Bar geometry: 100 µm long, 10 µm × 10 µm cross-section, sitting on a heat sink.
  const BAR = { L: 0.01, A: 1e-6, Rth: 60, tauTh: 1.2 };   // cm, cm², K/W, s

  // Silicon band gap (Varshni) and intrinsic density
  const Eg = T => 1.170 - 4.73e-4 * T * T / (T + 636);
  const Nc = T => 2.86e19 * Math.pow(T / 300, 1.5);
  const Nv = T => 3.10e19 * Math.pow(T / 300, 1.5);
  const ni = T => Math.sqrt(Nc(T) * Nv(T)) * Math.exp(-Eg(T) / (2 * kB * T));

  // Phosphorus donor level; it shrinks with doping and vanishes at the Mott transition
  // (≈3.7e18 cm⁻³), where the donors merge into a band and stop freezing out.
  const Ed = N => Math.max(0, 0.045 * (1 - Math.cbrt(N / 3.7e18)));

  // Charge neutrality n = p + N_D⁺ with partial donor ionisation (freeze-out at low T).
  function carriers(N, T) {
    const ni2 = ni(T) ** 2, nc = Nc(T), ed = Ed(N);
    const Nplus = n => N / (1 + 2 * (n / nc) * Math.exp(ed / (kB * T)));
    const f = n => n - ni2 / n - Nplus(n);
    let lo = 0, hi = 22;                       // log10 n
    for (let i = 0; i < 80; i++) { const m = (lo + hi) / 2; if (f(10 ** m) > 0) hi = m; else lo = m; }
    const n = 10 ** ((lo + hi) / 2);
    return { n, p: ni2 / n, ion: N > 0 ? Nplus(n) / N : 1, ni: Math.sqrt(ni2) };
  }

  // Low-field mobility, Arora et al. (1982): lattice (phonon) scattering plus
  // ionised-impurity scattering, as a function of doping and temperature. cm²/Vs.
  function mobility(N, T, hole) {
    const t = T / 300, a = 0.88 * Math.pow(t, -0.146);
    if (!hole) {
      const mn = 88 * Math.pow(t, -0.57), ml = 7.4e8 * Math.pow(T, -2.33);
      return { mu: mn + ml / (1 + Math.pow(N / (1.26e17 * Math.pow(t, 2.4)), a)), lat: mn + ml };
    }
    const mn = 54.3 * Math.pow(t, -0.57), ml = 1.36e8 * Math.pow(T, -2.23);
    return { mu: mn + ml / (1 + Math.pow(N / (2.35e17 * Math.pow(t, 2.4)), a)), lat: mn + ml };
  }

  // Velocity saturation, Caughey–Thomas form with Canali's temperature fits. cm/s.
  const vsat = (T, hole) => (hole ? 0.84 : 1) * 2.4e7 / (1 + 0.8 * Math.exp(T / 600));
  // (Canali's β fits only cover 300–430 K, so hold them at their 250 K value below that.)
  const beta = (T, hole) => { const t = Math.max(250, T) / 300; return hole ? 1.213 * Math.pow(t, 0.17) : 1.109 * Math.pow(t, 0.66); };
  function drift(mu, E, T, hole) {
    const vs = vsat(T, hole), b = beta(T, hole), x = mu * E / vs;
    return mu * E / Math.pow(1 + Math.pow(x, b), 1 / b);
  }

  // Everything about the bar at temperature T and applied voltage V.
  function state(N, T, V) {
    const c = carriers(N, T), me = mobility(N, T, false), mh = mobility(N, T, true);
    const E = Math.abs(V) / BAR.L;                         // V/cm
    const vn = drift(me.mu, E, T, false), vp = drift(mh.mu, E, T, true);
    const J = q * (c.n * vn + c.p * vp);                   // A/cm²
    const I = J * BAR.A;
    const sigma0 = q * (c.n * me.mu + c.p * mh.mu);        // low-field conductivity, S/cm
    const R0 = BAR.L / (sigma0 * BAR.A);
    return { ...c, N, T, V, E, muN: me.mu, muNlat: me.lat, muP: mh.mu, muPlat: mh.lat, vn, vp, I,
      R: I > 0 ? Math.abs(V) / I : R0, R0, P: Math.abs(V) * I, vsat: vsat(T, false) };
  }

  // Low-field resistance at temperature T (no self-heating).
  const R0 = (N, T) => state(N, T, 1e-3).R0;

  // Joule heating: the bar relaxes toward T_sink + P·R_th.
  function makeBar(N, Tsink, V) { return { N, Tsink, V, T: Tsink, dTdt: 0, burnt: false }; }
  function stepBar(b, dt) {
    let s = state(b.N, b.T, b.V);
    const Teq = b.Tsink + s.P * BAR.Rth;
    const k = 1 - Math.exp(-dt / BAR.tauTh);
    const dT = (Teq - b.T) * k;
    b.T += dT; b.dTdt = dt > 0 ? dT / dt : 0;
    if (b.T > 750) b.burnt = true;
    return state(b.N, b.T, b.V);
  }
  // Steady state, found by iterating the thermal balance (for checking levels).
  function steady(N, Tsink, V) {
    let T = Tsink, s;
    for (let i = 0; i < 400; i++) { s = state(N, T, V); const Tn = Tsink + s.P * BAR.Rth; T += 0.3 * (Tn - T); if (T > 2000) break; }
    return state(N, T, V);
  }

  // ---------- microscopic picture for the animation ----------
  // Per-carrier scattering time from the (field-dependent) mobility: τ = μ m*/q.
  // Positions in nm, time in ps, velocities in nm/ps.
  function micro(s) {
    const mk = hole => {
      const m = (hole ? MH : ME) * m0;
      const mu = hole ? s.muP : s.muN, mulat = hole ? s.muPlat : s.muNlat;
      const v = hole ? s.vp : s.vn;                        // cm/s
      const E = s.E * 100;                                 // V/m
      const muEff = s.E > 0 ? v / s.E : mu;                // cm²/Vs
      const tau = muEff * 1e-4 * m / q * 1e12;             // ps
      const acc = q * E / m * 1e-15;                       // nm/ps²
      // carrier heating: energy relaxation time ~0.3 ps
      const Te = Math.min(4000, s.T + (2 / 3) * q * (muEff * 1e-4) * E * E * 0.3e-12 / kJ);
      const vth = Math.sqrt(kJ * Te / m) * 1e-3;           // per-axis thermal speed, nm/ps
      // share of scatterings by mechanism (Matthiessen's rule)
      const rl = 1 / mulat, ri = Math.max(0, 1 / mu - 1 / mulat), rf = Math.max(0, 1 / muEff - 1 / mu);
      const tot = rl + ri + rf;
      return { tau, acc, vth, Te, drift: v * 1e-5,          // cm/s → nm/ps
        pLat: rl / tot, pImp: ri / tot, pHot: rf / tot, mfp: vth * Math.sqrt(2) * tau };
    };
    return { e: mk(false), h: mk(true) };
  }

  const api = { BAR, Eg, ni, Ed, carriers, mobility, vsat, drift, state, R0, makeBar, stepBar, steady, micro };
  if (typeof module !== 'undefined') module.exports = api; else root.Carrier = api;
})(this);
