# NeutronStars

End-to-end spectral image cube pipeline for neutron stars, built from first-principles physics in Julia.

## Goal

Compute `I(x, y, ν)` — specific intensity as a function of image-plane position and frequency — from the Schrödinger equation to the final pixel. Every equation traces to a locally-stored published paper. Every approximation is named and recorded.

## Pipeline

1. **Neutron star structure** — TOV solver with BSk EOS analytical fits (Potekhin+ 2013)
2. **Opacities** — Free-free (magnetic and non-magnetic), Thomson scattering, dielectric tensor normal modes (Potekhin & Chabrier 2003)
3. **Atmosphere RT** — Feautrier solver with Rybicki temperature correction, two-mode magnetic extension (Haakonsen+ 2012, Suleimanov+ 2009)
4. **Surface model** — Dipole magnetic field geometry, Greenstein-Hartke temperature map
5. **GR ray tracing** — Schwarzschild geodesics via elliptic integrals (Pechenick+ 1983, Beloborodov 2002)
6. **Rendering** — CIE 1931 colorimetry, sRGB, Reinhard tone mapping

## RX J1856.5−3754 — The Nearest Neutron Star

Rendered with real radiative transfer atmosphere spectra. M = 1.4 M☉, R = 14 km, B = 10¹³ G, T_pole = 7×10⁵ K, distance = 123 pc.

### X-ray rotation (one full period)

The hot magnetic pole sweeps into and out of view as the star rotates, producing the observed X-ray pulsation. False colour: red = 0.1–0.5 keV, green = 0.5–2 keV, blue = 2–10 keV.

![Rotation animation](output/rxj1856_rotation.gif)

### Spectral sweep (far-IR through hard X-ray)

Each frame images the star at a single photon energy, sweeping from 0.001 keV (far infrared) to 2.5 keV (hard X-ray). The hot pole is invisible at low energies but dominates in soft X-rays — explaining why the optical counterpart (HST, V≈25.7) shows weak pulsations while the X-ray lightcurve (Chandra/XMM) pulses strongly.

![Spectral sweep animation](output/rxj1856_spectral_sweep.gif)

## Current Status

- **Phase 2** (tracer bullet): Complete — TOV, ray tracer, surface model, colorimetry, spectral image cube with modified blackbody
- **Phase 3a** (non-magnetic atmosphere): Complete — temperature profile matches McPHAC within 1.2%, flux conservation F/σT⁴ = 0.99, frequency-adaptive depth grid (2× faster)
- **Phase 3b** (magnetic atmosphere): Working — B=10¹² converges (F/σT⁴=1.03), B=10¹⁴ converges (F/σT⁴=1.01), spectral hardening and proton cyclotron features verified against Suleimanov+(2009)
- **Phase 3c** (SpectralImageCube v2): Complete — real atmosphere spectra replace modified blackbody, 512×512×50 render in 4.6s

## Quick Start

```bash
julia --project=. -e 'using NeutronStar; println("OK")'

# Run non-magnetic atmosphere solver
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
result = solve_atmosphere(1e6, 2e14, gaunt; K=50, M=8, N=200, verbose=true)'
```

## Physics Report

A 30-page [physics report](docs/physics_report.pdf) documents every equation in the pipeline, cross-referenced against the source publications with equation numbers. It also serves as a code review — three discrepancies found (all minor: a misleading comment, an older solar mass constant, and the Gaunt factor table dependency).

## References

- Haakonsen et al. (2012) ApJ 749:52 — McPHAC (verification target)
- Potekhin & Chabrier (2003) ApJ 585:955 — Magnetic free-free opacities
- Suleimanov, Potekhin & Werner (2009) A&A 500:891 — Magnetic atmosphere scheme
- Pechenick, Ftaclas & Cohen (1983) ApJ 274:846 — NS image in Schwarzschild
- Beloborodov (2002) ApJ 566:L85 — Simplified geodesic formulae

## License

GPL-3.0. See [LICENSE](LICENSE).
