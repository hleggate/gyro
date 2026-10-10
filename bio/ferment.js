// Stirred-tank fermentation of Lactobacillus casei on glucose.
//
// Amounts are tracked as totals in the vessel (g or mol) so that feeding and
// base addition dilute everything correctly; concentrations are amount / V.
// Time is in hours.
//
// Growth (Monod with substrate and product inhibition, Baranyi lag):
//   mu = mu_opt * g(T) * g(pH) * S/(Ks + S + S^2/Ki) * (1 - P/Pc)^n * q/(1+q)
// g(T) is Rosso's cardinal temperature model, g(pH) the cardinal pH model.
// Lactic acid (Luedeking-Piret): dP/dt = (alpha * mu + beta) * Xv
// Death: acid stress below pH ~4.4, starvation once the glucose is gone,
// and heat above ~44 C. Dead cells still scatter light, so OD keeps them.
// pH: charge balance of lactic acid, the medium's buffers (acetate,
// phosphate, peptides) and the NaOH added by the pH-stat, solved by bisection.
(function (root) {
  const MU_OPT = 0.75;          // 1/h at 37-38 C, pH 6.3
  const KS = 0.4, KI = 160;     // g/L
  const PC = 100, NP = 0.6;     // lactate that stops growth (g/L)
  const PCP = 125;              // ... and production
  const ALPHA = 4.5, BETA = 0.15;  // g lactate per g cells; g/(g h)
  const YXS = 0.8;              // g cells per g glucose (for the biomass part)
  const YPS = 0.95;             // g lactate per g glucose
  const CFU_PER_G = 2.5e9;      // CFU/mL per g/L dry cells
  const OD_PER_G = 1 / 0.35;
  const MW = 90.08;             // lactic acid

  const T_MIN = 12, T_OPT = 38, T_MAX = 46;
  const PH_MIN = 4.0, PH_OPT = 6.4, PH_MAX = 8.6;

  // Rosso cardinal temperature model with inflection (CTMI)
  function gT(T) {
    if (T <= T_MIN || T >= T_MAX) return 0;
    const a = (T - T_MAX) * (T - T_MIN) ** 2;
    const b = (T_OPT - T_MIN) * ((T_OPT - T_MIN) * (T - T_OPT) - (T_OPT - T_MAX) * (T_OPT + T_MIN - 2 * T));
    return Math.max(0, a / b);
  }
  // cardinal pH model
  function gPH(pH) {
    if (pH <= PH_MIN || pH >= PH_MAX) return 0;
    const a = (pH - PH_MIN) * (pH - PH_MAX);
    return Math.max(0, a / (a - (pH - PH_OPT) ** 2));
  }

  // weak acids in the medium, mol/L at the start: [conc, pKa]
  const BUFFERS = [[0.03, 4.76], [0.015, 7.2], [0.03, 6.3], [0.02, 4.2]];
  const PKA_LA = 3.86, KW = 1e-14;
  // net charge of everything except the strong ions, at a given pH
  function anions(pH, c, lac) {
    const H = 10 ** -pH;
    let q = KW / H - H;
    for (const [ci, pk] of BUFFERS) { const K = 10 ** -pk; q += c * ci * K / (K + H); }
    const K = 10 ** -PKA_LA;
    return q + lac * K / (K + H);
  }
  // strong ion difference (mol/L) that gives pH0 in fresh medium
  const sid0 = pH0 => anions(pH0, 1, 0);

  function solvePH(f) {
    const c = f.V0 / f.V;                 // dilution of the original medium
    const lac = f.P / MW / f.V;           // mol/L total lactic acid
    const sid = (f.sid * f.V0 + f.base) / f.V;
    let lo = 2, hi = 11;
    for (let i = 0; i < 40; i++) { const m = (lo + hi) / 2; if (anions(m, c, lac) < sid) lo = m; else hi = m; }
    return (lo + hi) / 2;
  }

  // V0 (L), S0 (g/L glucose), pH0
  function makeReactor(o = {}) {
    const V0 = o.V0 || 1.5, pH0 = o.pH0 || 6.5;
    const f = {
      V0, V: V0, Vmax: o.Vmax || 2.0, X: 0, Xv: 0, S: (o.S0 ?? 20) * V0, P: 0, base: 0,
      sid: sid0(pH0), q: 0, t: 0, T: o.T ?? 37, setT: o.T ?? 37, pHset: 6.3,
      baseOn: false, feedOn: false, baseFlow: 0, inoculated: false, full: false, fed: 0,
      baseRate: 0.45,       // mol/h max (10 M NaOH, 45 mL/h)
      feedRate: 0.04,       // L/h of 600 g/L glucose
      feedConc: 600,
      peakCFU: 0,
    };
    f.pH = solvePH(f);
    return f;
  }

  function inoculate(f, gPerL = 0.02, lag = 1.5) {
    if (f.inoculated) return;
    f.X = f.Xv = gPerL * f.V;
    f.q = 1 / (Math.exp(lag * MU_OPT) - 1);
    f.inoculated = true;
  }

  function rates(f) {
    const S = f.S / f.V, P = f.P / f.V, T = f.T, pH = f.pH;
    const fS = S / (KS + S + S * S / KI);
    const fP = Math.pow(Math.max(0, 1 - P / PC), NP);
    const lagF = f.q / (1 + f.q);
    const mu = MU_OPT * gT(T) * gPH(pH) * fS * fP * lagF;
    const prod = gT(T) * Math.sqrt(gPH(pH)) * fS * Math.max(0, 1 - P / PCP);
    const qP = ALPHA * mu + BETA * prod;
    let kd = 0.004;
    if (pH < 4.4) kd += 0.7 * (4.4 - pH);
    if (S < 0.05) kd += 0.05 * (1 - S / 0.05);
    if (T > 44) kd += 0.6 * (T - 44);
    if (P > 90) kd += 0.01 * (P - 90);
    return { mu, qP, kd, fS };
  }

  function step(f, dt) {
    const n = Math.max(1, Math.ceil(dt / 0.01)), h = dt / n;
    for (let s = 0; s < n; s++) {
      f.T += (f.setT - f.T) * Math.min(1, h / 0.25);   // water jacket, ~15 min
      f.full = f.V >= f.Vmax - 1e-6;
      // pumps
      if (f.feedOn && !f.full) {
        const dv = f.feedRate * h; f.V += dv; f.S += f.feedConc * dv; f.fed += f.feedConc * dv;
      }
      f.baseFlow = 0;
      if (f.baseOn && !f.full && f.pH < f.pHset) {
        f.baseFlow = f.baseRate * Math.min(1, (f.pHset - f.pH) / 0.05);
        const dn = f.baseFlow * h;
        f.base += dn; f.V += dn / 10;
      }
      if (f.inoculated) {
        const r = rates(f);
        const g = r.mu * f.Xv * h, dP = r.qP * f.Xv * h;
        let dS = dP / YPS + g / YXS;
        const scale = dS > f.S ? f.S / dS : 1;           // can't eat glucose that isn't there
        f.Xv += g * scale - r.kd * f.Xv * h;
        f.X += g * scale;
        f.P += dP * scale; f.S = Math.max(0, f.S - dS * scale);
        f.q += MU_OPT * gT(f.T) * f.q * h;
      }
      f.pH = solvePH(f);
      f.t += h;
    }
    const c = cfu(f); if (c > f.peakCFU) f.peakCFU = c;
  }

  const cfu = f => f.Xv / f.V * CFU_PER_G;
  const od = f => f.X / f.V * OD_PER_G;
  const conc = f => ({ S: f.S / f.V, P: f.P / f.V, X: f.X / f.V, Xv: f.Xv / f.V });

  root.Ferment = { makeReactor, inoculate, step, rates, cfu, od, conc, gT, gPH, solvePH };
  if (typeof module !== 'undefined') module.exports = root.Ferment;
})(typeof window !== 'undefined' ? window : globalThis);
