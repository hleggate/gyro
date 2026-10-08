// Charged particles in a magnetic field pointing out of the screen (+z).
// Boris pusher, periodic boundaries. Units: pixels and simulation steps.
(function (root) {
  const SPECIES = {
    ion:      { qm:  1, speed: 3, color: '#ff9f43' },
    electron: { qm: -4, speed: 6, color: '#4fc3f7' },
  };

  function makeSim(W, H) {
    return {
      W, H,
      B0: 0.05,          // so ion gyrofrequency = 0.05 rad/step
      gradB: 0,          // fractional change of B across the screen height
      Ex: 0, Ey: 0,      // electric field (set by tilt or drag)
      particles: [],
    };
  }

  function Bat(sim, y) {
    // stronger towards the top of the screen (small y)
    const b = sim.B0 * (1 + sim.gradB * (0.5 - y / sim.H));
    return Math.max(b, sim.B0 * 0.05);
  }

  function addParticle(sim, kind, x, y) {
    const s = SPECIES[kind];
    const a = Math.random() * Math.PI * 2;
    const v = s.speed * (0.7 + 0.6 * Math.random());
    sim.particles.push({ kind, qm: s.qm, color: s.color, x, y,
      vx: v * Math.cos(a), vy: v * Math.sin(a), px: x, py: y });
  }

  function step(sim, dt) {
    const { W, H, Ex, Ey } = sim;
    for (const p of sim.particles) {
      p.px = p.x; p.py = p.y;
      const h = 0.5 * p.qm * dt;
      // half electric kick
      let vx = p.vx + h * Ex, vy = p.vy + h * Ey;
      // magnetic rotation
      const t = h * Bat(sim, p.y);
      const s = 2 * t / (1 + t * t);
      const v1x = vx + vy * t, v1y = vy - vx * t;
      vx += v1y * s; vy -= v1x * s;
      // second half electric kick
      p.vx = vx + h * Ex; p.vy = vy + h * Ey;
      p.x += p.vx * dt; p.y += p.vy * dt;
      p.wrapped = false;
      if (p.x < 0) { p.x += W; p.wrapped = true; } else if (p.x >= W) { p.x -= W; p.wrapped = true; }
      if (p.y < 0) { p.y += H; p.wrapped = true; } else if (p.y >= H) { p.y -= H; p.wrapped = true; }
    }
  }

  root.Gyro = { SPECIES, makeSim, addParticle, step, Bat };
  if (typeof module !== 'undefined') module.exports = root.Gyro;
})(typeof window !== 'undefined' ? window : globalThis);
