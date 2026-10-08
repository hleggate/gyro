# Gyro: pocket plasma

Ions and electrons gyrating in a magnetic field that points out of the screen. Tilt your phone (or drag a finger) to add an electric field and watch the E×B drift. Toggle ∇B to make the field stronger at the top, and the two species drift apart in opposite directions.

It's a static web app (HTML + canvas, no build step) served by GitHub Pages. On an iPhone, open it in Safari and use **Share → Add to Home Screen** to get it full-screen and offline.

- `physics.js`: Boris particle pusher, periodic boundaries
- `index.html`: drawing, touch and tilt input
- `sw.js`: offline cache (bump `VERSION` when files change)

## Laser Lab (`laser/`)

A second game: cut micron-sized channels in steel with a pulsed laser, seen in cross-section. Each pulse removes depth δ·ln(F/F_th) wherever the Gaussian beam is above the ablation threshold, reduced by cos θ on sloped walls, which is why real laser-cut walls taper. Femtosecond pulses cut cleanly; nanosecond pulses melt and leave burrs and a heat-affected zone. Four levels plus free play.

- `laser/laser.js`: ablation model and scoring
- `laser/index.html`: game, levels and drawing
