# Gyro: pocket plasma

Ions and electrons gyrating in a magnetic field that points out of the screen. Tilt your phone (or drag a finger) to add an electric field and watch the E×B drift. Toggle ∇B to make the field stronger at the top, and the two species drift apart in opposite directions.

It's a static web app (HTML + canvas, no build step) served by GitHub Pages. On an iPhone, open it in Safari and use **Share → Add to Home Screen** to get it full-screen and offline.

- `physics.js`: Boris particle pusher, periodic boundaries
- `index.html`: drawing, touch and tilt input
- `sw.js`: offline cache (bump `VERSION` when files change)

## Laser Lab (`laser/`)

A second game: cut micron-sized channels in steel with a pulsed laser, seen in cross-section. Each pulse removes depth δ·ln(F/F_th) wherever the Gaussian beam is above the ablation threshold, reduced by cos θ on sloped walls, which is why real laser-cut walls taper. Femtosecond pulses cut cleanly; nanosecond pulses melt and leave burrs and a heat-affected zone. Four levels plus free play.

Beam controls, available on every level:

- **Beam: Gaussian / Top hat.** A beam shaper redistributes the same pulse energy into a flat-topped (order-16 super-Gaussian) profile. The peak drops to 0.58× the Gaussian's, so you need more energy to reach threshold, but the floor comes out flat and the walls steep. A small plot above the surface shows the fluence profile against the ablation threshold.
- **GDD.** Group delay dispersion chirps the femtosecond pulse: τ = τ₀√(1 + (4 ln2·GDD/τ₀²)²), with τ₀ = 300 fs. Past ~10 ps (the electron–phonon time in steel) heat spreads during the pulse, so the threshold rises and the metal starts to melt, sliding from femtosecond-like to nanosecond-like behaviour (burrs and a heat-affected zone). It does nothing to the narrowband nanosecond laser.

- `laser/laser.js`: ablation model and scoring
- `laser/index.html`: game, levels and drawing

## Dot Foundry (`qdots/`)

Grow CdSe quantum dots by hot injection and fill orders for exact emission colours. Nucleation follows classical theory (J ∝ exp(−B/ln²S)) in a short burst after each injection; dots grow by diffusion-limited growth with Gibbs–Thomson solubility, so small dots dissolve and big ones grow. Emission wavelength comes from the Yu et al. (2003) CdSe sizing curve plus a Stokes shift. Temperature sets how many nuclei form (and so the final size); a second hot injection nucleates a new population, while slow drip feeding grows the existing dots. Scored on colour purity: the fraction of emission within ±15 nm of the order.

- `qdots/qd.js`: nucleation, growth and spectrum model
- `qdots/index.html`: flask, spectrum and orders

## Lacto Works (`bio/`)

Run a 2 L stirred-tank bioreactor growing *Lactobacillus casei* on glucose, and fill orders for live cultures and lactic acid. Growth is Monod kinetics with a Baranyi lag phase, slowed by substrate and lactic acid inhibition, a cardinal temperature model (optimum 38 °C, death above ~44 °C) and a cardinal pH model. Lactic acid comes from growth plus a non-growth term (Luedeking–Piret). The pH is solved from a charge balance of lactic acid, the medium's buffers and the NaOH added by the pH-stat; a bromocresol purple tint shows it in the broth. Cells die in acid below pH ~4.4 and when the glucose runs out, so harvest timing matters. The feed (600 g/L glucose) and base both add volume, and the vessel only holds 2 L. Four orders (starter culture, probiotic capsules, high-density fed-batch, a lactic acid titre for PLA) plus free play.

- `bio/ferment.js`: growth, acid production, death and pH model
- `bio/index.html`: reactor, growth chart and orders

## Free Path (`electrons/`)

Push electrons through a phosphorus-doped silicon bar, 100 µm long and 10 µm across. A 400 nm window shows a Monte Carlo of the carriers: free flights accelerated by the field, interrupted by scattering off phonons (white), ionised donors (orange) and, in strong fields, optical-phonon emission (purple). The scattering time comes from the mobility (τ = μm*/q), so the drift you see matches the current. Low-field mobility follows Arora et al. (1982) in doping and temperature, velocity saturation the Caughey–Thomas form with Canali's fits, and the carrier densities a charge balance with a 45 meV donor level (so donors freeze out when cold) plus intrinsic electron–hole pairs (so lightly doped silicon conducts more when hot). The bar heats itself through a 60 K/W heat sink and burns out above 750 K. Five levels (hit a current, make a resistor, beat the speed limit, work at 77 K, make a temperature-proof resistor) plus free play.

- `electrons/carrier.js`: carrier density, mobility, saturation and heating model, plus the Monte Carlo parameters
- `electrons/index.html`: crystal view, charts and levels
