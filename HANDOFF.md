# Session Handoff — 2026-03-24

## Project Overview

Building a physically exact neutron star visualisation pipeline in Julia. Master goal (from `NEUTRON_STAR_CLAUDE_MD_v2.md`): spectral image cube I(x,y,ν) from first principles. Every equation traceable to a locally-stored published paper.

GitHub: https://github.com/tobiasosborne/neutronstars (public, GPL-3.0)

---

## What's Working

### Phase 2 (Tracer Bullet) — COMPLETE
- TOV solver with BSk EOS (verified against Potekhin+ 2013)
- Schwarzschild ray tracer (verified against Pechenick+ 1983, Beloborodov 2002)
- Dipole surface model (Greenstein-Hartke temperature map)
- CIE 1931 colorimetry + sRGB rendering (verified: 5778K solar white)
- Full 512×512 renders working

### Phase 3a (Non-Magnetic Atmosphere) — COMPLETE
- Feautrier RT solver with Rybicki (1971) temperature correction
- **Temperature profile matches McPHAC within 1.2%** at T_eff=10⁶ K
- **Flux conservation F/σT⁴ = 0.99** (verified)
- **Multi-parameter verification**: 10⁶ K (1.2%), 10⁵·³ K (5.4%), 10⁶·⁵ K (0.96%)
- **Frequency-adaptive depth grid**: `solve_feautrier_all_adaptive()` — 2× faster, <0.3% difference

### Phase 3b (Magnetic Atmosphere) — WORKING
- `magnetic_atmosphere.jl`: two-mode Feautrier + Rybicki correction + uniform flux correction
- `mode_absorption()` and `mode_scattering()` properly separated in `magnetic_modes.jl`
- **B=0 limit**: matches non-magnetic within 0.06%
- **B=10¹² G**: converges in 60 iterations, F/σT⁴ = 1.026
- **B=10¹⁴ G**: converges in ~150 iterations, F/σT⁴ = 1.014
- Spectral hardening: 10× blackbody at 1 keV (B=10¹²), proton cyclotron feature at 0.635 keV (B=10¹⁴)
- Qualitative match to SPW2009 (quantitative differences from missing vacuum polarization)

### Phase 3c (SpectralImageCube v2) — COMPLETE
- `AtmosphereGrid` module: pre-computes atmosphere spectra at grid of (T_eff, B) values
- `SpectralImageCube` struct matching the master spec in `NEUTRON_STAR_CLAUDE_MD_v2.md`
- `render_spectral_cube()` replaces modified-blackbody placeholder with real RT spectra
- 256×256 render with 50 frequencies in 1.3s (grid build: 14s for B=0, ~12min with B=10¹²)
- Output: `output/ns_v2_256_{true,xray}.ppm`

### Magnetic Opacity Infrastructure — COMPLETE
- Dielectric tensor (Stix parameters S,D,P), q-based polarization vectors
- Magnetic free-free cross-sections σ_α for α = -1, 0, +1 (P&C 2003)
- Magnetic Coulomb logarithm with Landau quantization (Eq. 44)
- Mode decomposition: κ_j(ν, θ_B, B, T, ρ) for j=X,O
- `mode_absorption()` + `mode_scattering()` properly separated
- Angle-averaged Rosseland means K_∥, K_⊥ (verified: 0.17 dex at low density)

---

## What Needs Doing Next (Priority Order)

### 1. Improve Magnetic Atmosphere Convergence
The uniform flux correction works but is slow (60 iterations for B=10¹²). SPW2009 uses a three-part correction (Λ-iteration + Avrett-Krook + surface), which would converge faster and produce better T profiles. The T profile at B=10¹⁴ is flatter than physical (uniform correction limitation).

### 2. Wire Adaptive Feautrier into RT Iteration
`solve_feautrier_all_adaptive()` exists but isn't yet used in `solve_atmosphere()` or `solve_magnetic_atmosphere()`. Replacing `solve_feautrier_all` with the adaptive version in the iteration loop should improve spectral accuracy from ~7% to <1%.

### 3. Add Vacuum Polarization
SPW2009 shows vacuum polarization significantly affects the proton cyclotron absorption feature. Implement mode conversion probability P_jump (SPW2009 Eq. 16-17) using van Adelsberg & Lai (2006).

### 4. HDF5 Atmosphere Grid Storage
Currently the atmosphere grid is computed fresh each session. Store as HDF5 with provenance metadata for reuse. Requires adding HDF5.jl dependency.

### 5. (Phase 4) Extensions
In order: OCP Monte Carlo, partial ionisation, vacuum QED mode conversion, **Kerr ray tracing** (see Lyr.jl below), magnetosphere volume rendering. See `NEUTRON_STAR_CLAUDE_MD_v2.md` Phase 4.

### Lyr.jl — Kerr Ray Tracer for Phase 4
`/home/tobiasosborne/Projects/Lyr.jl/` has a production-ready Kerr ray tracer:
- Full Boyer-Lindquist metric with analytic Christoffels (10-20× faster than ForwardDiff)
- 169 tests passing, ISCO formulas, frame dragging
- Modular GR module: `src/GR/metrics/kerr.jl` + `src/GR/integrator.jl` + `src/GR/camera.jl`
- Can be extracted without the volume rendering / OpenVDB dependencies

---

## Architecture

### Module Dependency Order (NeutronStar.jl)
```
PhysicalConstants → BSkEOS → TOVSolver → DipoleModel
                  → GauntFactor → HydrogenOpacity → MagneticCoulomb → MagneticFF
                  → DielectricTensor → BlackbodyAtmosphere → MagneticModes
                  → AtmosphereStructure → FeautrierSolver → TemperatureCorrection
                  → RTAtmosphere → MagneticAtmosphere
                  → SchwarzschildTracer → CIE_sRGB → AtmosphereGrid → Renderer
```

### Key Files and What They Do

**Atmosphere:**
| File | Purpose | Status |
|------|---------|--------|
| `src/atmosphere/feautrier.jl` | Block-tridiagonal Feautrier RT solver + adaptive variant | Working, verified vs McPHAC |
| `src/atmosphere/temp_correction.jl` | Rybicki (1971) global temperature correction | Working, verified vs McPHAC |
| `src/atmosphere/rt_atmosphere.jl` | Non-magnetic atmosphere driver | Working |
| `src/atmosphere/magnetic_atmosphere.jl` | Two-mode magnetic atmosphere + flux correction | Working for B=0 to 10¹⁴ |
| `src/atmosphere/atm_structure.jl` | Column structure, Eddington T profile, opacity computation | Working |
| `src/atmosphere/blackbody.jl` | Planck function B_ν and modified blackbody | Working |

**Opacities:**
| File | Purpose | Status |
|------|---------|--------|
| `src/opacity/hydrogen_ff.jl` | B=0 free-free opacity (Haakonsen Eq. 12) + Thomson | Working |
| `src/opacity/gaunt_ff.jl` | Gaunt factor interpolation (Sutherland 1998 table) | Working |
| `src/opacity/magnetic_ff.jl` | Magnetic free-free σ_α cross-sections (P&C 2003) | Working |
| `src/opacity/magnetic_modes.jl` | Mode decomposition + `mode_absorption`/`mode_scattering` | Working |
| `src/opacity/dielectric_tensor.jl` | Cold-plasma Stix parameters, polarization weights | Working |
| `src/opacity/coulomb_magnetic.jl` | Magnetic Coulomb logarithm (Eq. 44) | Working |

**Pipeline:**
| File | Purpose | Status |
|------|---------|--------|
| `src/pipeline/atmosphere_grid.jl` | Pre-computed atmosphere spectrum grid + interpolation | NEW |
| `src/pipeline/render.jl` | SpectralImageCube + render_spectral_cube + RGB post-processing | Updated |

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
refs/code/McPHAC/OUT_T6_g2e14/    ← T_eff=10⁶, g=2×10¹⁴ (converged, N=400)
refs/code/McPHAC/OUT_T5.3_g1e14/  ← T_eff=10⁵·³, g=10¹⁴ (converged, N=400)
refs/code/McPHAC/OUT_T6.5_g2e14/  ← T_eff=10⁶·⁵, g=2×10¹⁴ (converged, N=400)
refs/code/McPHAC/gffgu.dat        ← Sutherland Gaunt factor table (used by both codes)
```

---

## Critical Physics Conventions

1. **Per-mode Planck source = B_ν/2**: Unpolarized thermal emission splits equally between X and O modes.
2. **Feautrier surface BC = pure radiation**: No local source (Q=0) and no scattering at surface (i=1).
3. **Rybicki surface = zero coupling**: `U_k[1]=0`, `K_k[1]=0`.
4. **EOS**: ρ = m_p P / (2 k_B T) for fully ionized hydrogen (μ=0.5).
5. **B̄ in Rybicki is GLOBAL**: Sum over both modes. Do NOT compute per-mode.
6. **Flux diagnostic**: F = 2π ∫₀¹ I(μ) μ dμ (not 4π).
7. **Uniform flux correction**: ΔT/T = -0.1×(flux_ratio - 1)/(4×flux_ratio). Applied when |F/σT⁴ - 1| > 2%.

---

## Quick Test Commands

```bash
# Compile check
julia --project=. -e 'using NeutronStar; println("OK")'

# Non-magnetic atmosphere (converges in ~10 iter, F/σT⁴ ≈ 1.0)
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
result = solve_atmosphere(1e6, 2e14, gaunt; K=50, M=8, N=200, max_iter=30, tol=1e-6, verbose=true)'

# Magnetic atmosphere B=10^12 (converges in ~60 iter, F/σT⁴ ≈ 1.02)
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
r = solve_magnetic_atmosphere(1e6, 2e14, 1e12, π/4, gaunt; K=50, M=8, N=100, tol=1e-4, max_iter=80, verbose=true)'

# SpectralImageCube v2 render (builds grid + renders 128×128)
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
grid = build_atmosphere_grid([5e5, 1e6, 2e6], [0.0], 2e14, gaunt; K=50, M=8, N=100, verbose=true)
params = NSParams(1.4, 12.0, 1e12, 1.5e6, 3e5, 0.3, π/3, 100.0, 1.0)
cube = render_spectral_cube(params, grid, 128; verbose=true)
tc, fc = render_cube_rgb(cube)
save_cube_ppm(tc, fc, 128, "output/ns_v2_128")'
```

---

## References

- `refs/haakonsen_2012_mcphac.pdf` — McPHAC paper (Appendix A = Rybicki method)
- `refs/potekhin_chabrier_2003_ff_opacity.pdf` — Magnetic opacities (Eqs. 25-52)
- `refs/suleimanov_potekhin_werner_2009_mag_atm.pdf` — Magnetic atmosphere models (SPW2009)
- `refs/potekhin_tables/hmagtab.txt` — Table format specification
- `NEUTRON_STAR_CLAUDE_MD_v2.md` — Master project specification (v2)
- **WARNING**: `refs/potekhin_chabrier_ho_2014_opacities.pdf` is the WRONG PAPER (Boev & Kovalev 2014 about exciton BEC). Needs replacement.
