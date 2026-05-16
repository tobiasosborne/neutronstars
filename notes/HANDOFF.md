# Session Handoff — 2026-05-16

## Project Overview

Physically traceable neutron-star atmosphere + rendering pipeline in Julia.
Master goal (from `CLAUDE.md`): spectral image cube
`I(x, y, ν)` from first principles, every equation traceable to a locally
stored published paper. Current active work is on magnetic atmospheres and
opacities, validated against local papers and digitized published figures.

GitHub: https://github.com/tobiasosborne/neutronstars (public, GPL-3.0)

This handoff covers the merge of two parallel work streams:
- **Local (Phase 3b + 3c)**: SpectralImageCube v2, AtmosphereGrid, RX J1856
  renders, 30-page physics report.
- **Remote (Suleimanov Fig 2 validation branch)**: free-free opacity bug
  fixes (Eq. 44e), vacuum-polarized normal-mode vectors, magnetic
  atmosphere initial-guess + flux-correction overhaul, Suleimanov Fig 2
  digitisation.

---

## Current Status

### Working And Verified

- **Phase 2 (tracer bullet)** — TOV (BSk EOS), Schwarzschild ray tracer,
  dipole surface model, CIE 1931 → sRGB rendering. Full 512×512 renders.
- **Phase 3a (non-magnetic atmosphere)** — Feautrier + Rybicki (1971).
  T profile matches McPHAC within 1.2% at T_eff=10⁶ K; F/σT⁴ = 0.99;
  multi-parameter verification at 10⁵·³, 10⁶, 10⁶·⁵ K.
- **Phase 3b (magnetic atmosphere)** — two-mode Feautrier with separated
  absorption / scattering, magnetic-Eddington initial T profile, adaptive
  column depth to diffusion cutoff, damped ad-hoc grey-body flux scaling
  (not the real SPW09 Avrett-Krook scheme; see notes/decisions.md D10),
  thermalisation-depth guard before LTE bottom BC.
  - Smoke test (`B=1×10¹² G`, `θ_B=π/4`): converges in ~60 iterations.
  - B=0 limit recovers non-magnetic to <0.1%.
- **Phase 3c (SpectralImageCube v2)** — `AtmosphereGrid` precomputes
  atmosphere spectra over (T_eff, B); `render_spectral_cube` replaces the
  modified-blackbody placeholder with real RT spectra; 256×256×50 render
  in ~1.3 s after the grid is built.
- **Magnetic free-free opacity** — matches the clean Suleimanov 2009
  Fig. 2 free-free branches outside resonance-marker bands. Eq. 44e
  transcription bug fixed; Eq. 44a integral extended to ∞; Eq. 52 proton
  damping split.
- **Vacuum-polarised normal modes** — `polarization_weights_vacuum` and
  vacuum coefficient fits (Potekhin, Lai & Chabrier 2004 A7–A9)
  implemented and regression-tested. `polarization_weights_full` remains
  as a cold-plasma compatibility wrapper.

### Still Not Physically Complete

- **Vacuum resonance mode conversion** is not implemented (van Adelsberg
  & Lai 2006 P_jump). Affects the proton cyclotron feature shape.
- **Scalar / local scattering approximation** — the Feautrier extinction
  term uses polarization-weighted scattering as a local cross section;
  full inter-mode angular redistribution is deferred.
- **Suleimanov Fig 2 clean validation** still shows a major mismatch at
  `θ_B = 45°` (~2.7–3.0 dex median residual on the upper branches at
  `B=10¹⁴ G`). `θ=90°` is good; `θ=5°` is mixed (~0.4–0.5 dex).
- **Render pipeline** still uses simplified placeholder behaviour in
  some paths — do not treat generated images as final physical products.

---

## Main Results From The Merged Sessions

### Magnetic RT Architecture (final form)

Files: `src/opacity/magnetic_modes.jl`, `src/atmosphere/magnetic_atmosphere.jl`,
`src/NeutronStar.jl`.

- Split magnetic opacity into:
  - `mode_absorption(j, ν, θ_B, B, T, ρ)` — true absorption (ff + pp)
  - `mode_scattering(j, ν, θ_B, B, T, ρ)` — scattering extinction
  - `mode_opacity = mode_absorption + mode_scattering`
- Scattering albedo `ρ_alb = σ_scat / κ_total` carried through Feautrier.
- Magnetic grey starting profile (`_magnetic_eddington_temperature`) uses
  the directional Rosseland mean of `effective_opacity`.
- Initial column depth grows automatically until every (mode, ν) reaches
  τ ≥ 80, capped at `y_max = 1e5 g cm⁻²`.
- Damped ad-hoc grey-body flux scaling in the iteration loop (NOT the
  real SPW09 Avrett-Krook scheme — see notes/decisions.md D10):
  `ΔT *= 1 + flux_damping · (flux_ratio^(-1/4) − 1)`, clamped to
  [0.9, 1.1]. Knobs are `flux_tol` (default 1e-2) and `flux_damping`
  (default 0.5).
- Convergence requires both `max|ΔT/T| < tol` and
  `|F/σT⁴ − 1| < flux_tol`.

### Magnetic Free-Free Ground Truth Fix

Files: `src/opacity/coulomb_magnetic.jl`, `src/opacity/magnetic_ff.jl`,
`test/runtests.jl`.

Ground truth: Potekhin & Chabrier 2003 Eq. 44a–e and Eq. 52.

- Eq. 44e correct form: `x_n = |u − n·βₑ| · √(1/4 + y/βₑ)`.
- Eq. 44a integral over `0 → ∞`; the previous `y=20` cutoff missed the
  long longitudinal-polarisation tail.
- Eq. 52 damping split now excludes proton-proton damping from the
  electron damping channel.
- Regression: low-energy Fig. 2 longitudinal point now has
  `Λ₀ ≈ 7.4148`. `verification/check_suleimanov_fig2_metrics.py` passes.

### Vacuum-Polarised Mode Vectors

Files: `src/opacity/dielectric_tensor.jl`, `src/opacity/magnetic_modes.jl`,
`test/runtests.jl`.

Ground truth: Potekhin, Lai & Chabrier 2004 vacuum coefficient fits;
van Adelsberg & Lai 2006 B-frame cyclic component formulae.

Functions added: `vacuum_coefficients`, `vacuum_polarization_parameter`,
`vacuum_mode_parameters`, `polarization_weights_vacuum`,
`polarization_weights_cold` (factored out of the old
`polarization_weights_full`, which is preserved as a compatibility
wrapper).

Regression target at `B=1×10¹⁴ G`, `E=1 keV`, `ρ=1 g cm⁻³`, `θ_B=45°`:
```
a_hat = -2.0563758036448108e-4
q_vac =  1.0012982615339183e-3
m_vac = -2.425084838547189e-4
beta  = 19.411715027031708
K1    = -0.025737455521326447
K2    = 38.84916750958474
```

### Suleimanov Fig 2 Digitisation And Validation

Files (verification):
- `verification/digitize_suleimanov_fig2.py`
- `verification/compute_suleimanov_fig2_opacities.jl`
- `verification/plot_suleimanov_fig2_opacity_comparison.py`
- `verification/plot_suleimanov_fig2_clean_validation.py`
- `verification/check_suleimanov_fig2_metrics.py`

Data + figures live under `verification/data/suleimanov_2009_fig2/` and
`verification/figures/suleimanov_2009_fig2/`.

Status:
- Free-free branches at `θ=5°` and `θ=90°` pass outside masked resonance
  bands.
- Clean-branch validation: `θ=90°` good, `θ=5°` mixed (~0.4–0.5 dex),
  `θ=45°` clearly bad (~2.7–3.0 dex on the upper branches).

**Do not claim Fig. 2 is fully matched.** The clean validation proves
the oblique-angle mismatch is real.

### SpectralImageCube v2 + AtmosphereGrid

Files: `src/pipeline/atmosphere_grid.jl`, `src/pipeline/render.jl`,
`src/NeutronStar.jl`.

- `AtmosphereGrid` module precomputes atmosphere spectra at a grid of
  `(T_eff, B)` and supports `lookup_spectrum`. Build cost: ~14 s for
  B=0, ~12 min with B=10¹² (60-iteration magnetic case).
- `SpectralImageCube` matches the master spec in
  `CLAUDE.md` (Section "Spectral image cube").
- `render_spectral_cube(params, grid, N)` replaces the modified-blackbody
  placeholder with real RT spectra; 256×256×50 in ~1.3 s.
- Outputs in `output/ns_v2_256_{true,xray}.ppm` and the RX J1856 GIFs
  (`output/rxj1856_rotation.gif`, `output/rxj1856_spectral_sweep.gif`).
- Frequency-adaptive Feautrier (`solve_feautrier_all_adaptive`) exists
  but is **not yet wired into the RT iteration loops** in
  `solve_atmosphere` / `solve_magnetic_atmosphere`.

### Physics Report

`docs/physics_report.{tex,pdf}` is a ~30-page equation-by-equation
review of the entire pipeline cross-referenced against the source
publications. It also flagged three minor discrepancies (a misleading
comment, an older solar mass constant, the Gaunt-factor table
dependency).

---

## Immediate Next Task

Focus on the `θ_B = 45°` clean Fig. 2 mismatch.

Recommended diagnostic:

1. Pick fixed energies away from masked bands, e.g.
   `0.03, 0.1, 0.3, 1.0, 3.0, 10.0 keV`.
2. For each energy at `T=7×10⁶ K`, `ρ=30 g cm⁻³`, `B=10¹⁴ G`,
   `θ_B=45°`, print:
   - `β`, `K_j`, `K_{z,j}`
   - `|e_-|²`, `|e_0|²`, `|e_+|²`
   - `sigma_ff_alpha`, `sigma_pp_alpha`, `sigma_scat_alpha`
   - final `mode_absorption`, `mode_scattering`
3. Compare directly against the paper equations, not against the
   plotted overlay.
4. Check whether the mismatch is caused by:
   - mode branch convention near oblique propagation,
   - B-frame cyclic weight sign/order,
   - missing scattering redistribution,
   - or the fact that the paper used van Adelsberg & Lai 2006 opacity
     formulae rather than our P&C component implementation.

### Secondary backlog (in priority order)

1. **Wire adaptive Feautrier into RT iteration** —
   `solve_feautrier_all_adaptive()` exists but isn't used by
   `solve_atmosphere` / `solve_magnetic_atmosphere`. Should improve
   spectral accuracy from ~7% to <1%.
2. **Vacuum resonance mode conversion** — P_jump (SPW2009 Eq. 16–17,
   van Adelsberg & Lai 2006).
3. **HDF5 atmosphere grid storage** — currently recomputed each session;
   adds an HDF5.jl dependency.
4. **Phase 4 extensions** — OCP Monte Carlo, partial ionisation, Kerr
   ray tracing (see Lyr.jl note below), magnetosphere volume rendering.

### Lyr.jl — Kerr Ray Tracer For Phase 4

`/home/tobiasosborne/Projects/Lyr.jl/` has a production-ready Kerr ray
tracer: full Boyer-Lindquist metric with analytic Christoffels (10–20×
faster than ForwardDiff), 169 tests passing, modular GR module
(`src/GR/metrics/kerr.jl`, `src/GR/integrator.jl`, `src/GR/camera.jl`).
Can be extracted without the volume-rendering / OpenVDB dependencies.

---

## Architecture

### Module Dependency Order (`src/NeutronStar.jl`)

```
PhysicalConstants → BSkEOS → TOVSolver → DipoleModel
                  → GauntFactor → HydrogenOpacity → MagneticCoulomb → MagneticFF
                  → DielectricTensor → BlackbodyAtmosphere → MagneticModes
                  → AtmosphereStructure → FeautrierSolver → TemperatureCorrection
                  → RTAtmosphere → MagneticAtmosphere
                  → SchwarzschildTracer → CIE_sRGB → AtmosphereGrid → Renderer
```

### Key Files

**Atmosphere:**
| File | Purpose | Status |
|------|---------|--------|
| `src/atmosphere/feautrier.jl` | Block-tridiagonal Feautrier RT solver + adaptive variant | Working, verified vs McPHAC |
| `src/atmosphere/temp_correction.jl` | Rybicki (1971) global temperature correction | Working, verified vs McPHAC |
| `src/atmosphere/rt_atmosphere.jl` | Non-magnetic atmosphere driver | Working |
| `src/atmosphere/magnetic_atmosphere.jl` | Two-mode magnetic atmosphere; magnetic-Eddington init; ad-hoc grey-body flux scaling (not real SPW09 Avrett-Krook, see decisions.md D10); thermalisation guard | Working for B=0 to 10¹⁴ G |
| `src/atmosphere/atm_structure.jl` | Column structure, Eddington T profile, opacity computation | Working |
| `src/atmosphere/blackbody.jl` | Planck function B_ν and modified blackbody | Working |

**Opacities:**
| File | Purpose | Status |
|------|---------|--------|
| `src/opacity/hydrogen_ff.jl` | B=0 free-free (Haakonsen Eq. 12) + Thomson | Working |
| `src/opacity/gaunt_ff.jl` | Gaunt factor interpolation (Sutherland 1998) | Working |
| `src/opacity/magnetic_ff.jl` | Magnetic free-free σ_α, proton-proton, scattering (P&C 2003) — Eq. 44e fixed | Working |
| `src/opacity/coulomb_magnetic.jl` | Magnetic Coulomb logarithm (Eq. 44, full ∞ integral) | Working |
| `src/opacity/magnetic_modes.jl` | `mode_absorption`, `mode_scattering`, `mode_opacity`, `effective_opacity` | Working |
| `src/opacity/dielectric_tensor.jl` | Stix parameters, cold + vacuum-polarised normal-mode weights | Working; θ_B=45° issue under investigation |

**Pipeline:**
| File | Purpose | Status |
|------|---------|--------|
| `src/pipeline/atmosphere_grid.jl` | Pre-computed atmosphere spectrum grid + interpolation | Working |
| `src/pipeline/render.jl` | `SpectralImageCube`, `render_spectral_cube`, RGB post-processing | Working |

**Stable (don't touch unless needed):**
| File | Purpose |
|------|---------|
| `src/constants.jl` | CODATA physical constants |
| `src/eos/bsk_eos.jl` | BSk analytical EOS fits |
| `src/eos/tov.jl` | TOV integrator |
| `src/surface/dipole.jl` | Dipole B-field + Greenstein-Hartke T map |
| `src/geodesics/schwarzschild.jl` | GR ray tracing via elliptic integrals |
| `src/colorimetry/cie_srgb.jl` | CIE 1931 → XYZ → sRGB pipeline |

### McPHAC Reference Data

```
refs/code/McPHAC/                  ← McPHAC source + binary (cloned from GitHub)
refs/code/McPHAC/OUT_T6_g2e14/     ← T_eff=10⁶,  g=2×10¹⁴ (converged, N=400)
refs/code/McPHAC/OUT_T5.3_g1e14/   ← T_eff=10⁵·³, g=10¹⁴ (converged, N=400)
refs/code/McPHAC/OUT_T6.5_g2e14/   ← T_eff=10⁶·⁵, g=2×10¹⁴ (converged, N=400)
refs/code/McPHAC/gffgu.dat         ← Sutherland Gaunt factor table (used by both codes)
```

`refs/code/McPHAC` is treated as a nested repo / submodule-style path
and shows persistent dirt in `git status`; do not commit changes into
it.

---

## Critical Physics Conventions

1. **Per-mode Planck source = B_ν/2** — unpolarised thermal emission
   splits equally between X and O modes.
2. **Feautrier surface BC = pure radiation** — no local source (Q=0)
   and no scattering at surface (i=1).
3. **Rybicki surface = zero coupling** — `U_k[1] = 0`, `K_k[1] = 0`.
4. **EOS** — `ρ = m_p · P / (2 · k_B · T)` for fully ionised hydrogen
   (μ = 0.5).
5. **B̄ in Rybicki is GLOBAL** — sum over both modes; do not compute
   per-mode.
6. **Flux diagnostic** — `F = 2π ∫₀¹ I(μ) μ dμ` (not 4π).
7. **Magnetic flux correction** — damped grey scaling
   `ΔT *= 1 + flux_damping · (flux_ratio^(-1/4) − 1)`, clamped to
   `[0.9, 1.1]`, applied only for `B > 0`.

---

## Commands To Reproduce Current Validation

```bash
# Compile check
julia --project=. -e 'using NeutronStar; println("OK")'

# Test suite
julia --project=. test/runtests.jl

# Non-magnetic atmosphere (converges in ~10 iter, F/σT⁴ ≈ 1.0)
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
result = solve_atmosphere(1e6, 2e14, gaunt; K=50, M=8, N=200, max_iter=30, tol=1e-6, verbose=true)'

# Magnetic atmosphere B=10¹² (converges in ~60 iter)
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
r = solve_magnetic_atmosphere(1e6, 2e14, 1e12, π/4, gaunt;
        K=50, M=8, N=100, max_iter=80, tol=1e-4,
        flux_tol=1e-2, flux_damping=0.5, verbose=true)'

# SpectralImageCube v2 render (builds grid + renders 128×128)
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
grid = build_atmosphere_grid([5e5, 1e6, 2e6], [0.0], 2e14, gaunt; K=50, M=8, N=100, verbose=true)
params = NSParams(1.4, 12.0, 1e12, 1.5e6, 3e5, 0.3, π/3, 100.0, 1.0)
cube = render_spectral_cube(params, grid, 128; verbose=true)
tc, fc = render_cube_rgb(cube)
save_cube_ppm(tc, fc, 128, "output/ns_v2_128")'

# Suleimanov Fig 2 validation
python3 verification/digitize_suleimanov_fig2.py
julia --project=. verification/compute_suleimanov_fig2_opacities.jl
python3 verification/plot_suleimanov_fig2_opacity_comparison.py
python3 verification/check_suleimanov_fig2_metrics.py
python3 verification/plot_suleimanov_fig2_clean_validation.py
```

---

## Important Caveats For The Next Agent

- `refs/ho_lai_2001_vacuum_pol.pdf` is **not** Ho & Lai — it is an
  unrelated quasar paper. Re-download from arXiv.
- `refs/potekhin_chabrier_ho_2014_opacities.pdf` is mislabeled — it is
  Boev & Kovalev 2014 (exciton BEC). Re-download from arXiv.
- arXiv source tarballs used during the Fig 2 session were downloaded
  to `/tmp` and are not committed:
  - `/tmp/van_adelsberg_lai_2006_src/ms14.tex`
  - `/tmp/potekhin_lai_chabrier_2004_src/kke.tex`
- `refs/code/McPHAC` appears dirty as a nested repo / submodule-style
  path — leave alone.
- Output PNG/PPM/GIF files under `output/` and LaTeX build artifacts
  under `docs/` (`*.aux`, `*.log`, `*.out`, `*.toc`) are working-tree
  noise that should probably move into `.gitignore`.

---

## References (Local PDFs)

- `refs/haakonsen_2012_mcphac.pdf` — McPHAC paper (Appendix A = Rybicki method)
- `refs/potekhin_chabrier_2003_ff_opacity.pdf` — Magnetic opacities (Eqs. 25–52)
- `refs/suleimanov_potekhin_werner_2009_mag_atm.pdf` — Magnetic atmosphere models (SPW2009)
- `refs/potekhin_tables/hmagtab.txt` — Table format specification
- `CLAUDE.md` — Master project specification
