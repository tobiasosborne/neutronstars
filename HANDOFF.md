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
- SpectralImageCube v1 with modified blackbody placeholder
- Full 512×512 renders working

### Phase 3a (Non-Magnetic Atmosphere) — WORKING
- Feautrier RT solver with Rybicki (1971) temperature correction
- **Temperature profile matches McPHAC within 1.2%** at T_eff=10⁶ K
- **Flux conservation F/σT⁴ = 0.99** (verified)
- **Multi-parameter verification**: 10⁶ K (1.2%), 10⁵·³ K (5.4%), 10⁶·⁵ K (0.96%)
- Spectrum matches within ~5-10% at peak (limited by fixed depth grid)

### Phase 3b (Magnetic Atmosphere) — FRAMEWORK ONLY
- `magnetic_atmosphere.jl`: two-mode Feautrier + Rybicki correction
- **B=0 limit verified**: matches non-magnetic solver within 0.06%
- **B>0 does NOT work yet** — see "Immediate Blocker" below

### Magnetic Opacity Infrastructure — COMPLETE
- Dielectric tensor (Stix parameters S,D,P), q-based polarization vectors
- Magnetic free-free cross-sections σ_α for α = -1, 0, +1 (P&C 2003)
- Magnetic Coulomb logarithm with Landau quantization (Eq. 44)
- Mode decomposition: κ_j(ν, θ_B, B, T, ρ) for j=X,O
- Angle-averaged Rosseland means K_∥, K_⊥ (verified: 0.17 dex at low density)

---

## Immediate Blocker: Magnetic Scattering Separation

The magnetic atmosphere solver (`magnetic_atmosphere.jl`) diverges for B>0 because `mode_opacity()` in `magnetic_modes.jl` returns the **total** opacity (absorption + scattering combined), but the Feautrier solver needs them **separated**:

- `κ_j` (absorption only) → enters thermal source term `S = (1-ρ)B + ρJ`
- `σ_j` (scattering only) → determines albedo `ρ_j = σ_j/(κ_j + σ_j)`
- `k_total_j = κ_j + σ_j` → determines optical depth scale

**Current `mode_opacity` combines them** via `σ_total_alpha` which sums `σ_ff_alpha + σ_pp_alpha + σ_scat_alpha`. The fix needs:

1. **Split `mode_opacity` into `mode_absorption` and `mode_scattering`** in `magnetic_modes.jl`:
   - `mode_absorption(j, ν, θ_B, B, T, ρ)` → uses only `σ_ff_alpha + σ_pp_alpha`
   - `mode_scattering(j, ν, θ_B, B, T, ρ)` → uses only `σ_scat_alpha`

2. **Update `_compute_magnetic_opacities!` in `magnetic_atmosphere.jl`** to call both and populate the `κ`, `σ`, `k_total`, `ρ_alb` arrays separately.

3. **Verify** that B→0 still recovers non-magnetic (both absorption and scattering should approach the B=0 values).

The underlying cross-section functions already exist in `magnetic_ff.jl`:
- `sigma_ff_alpha(α, ω, B, T, ρ)` — free-free absorption
- `sigma_pp_alpha(α, ω, B, T, ρ)` — proton-photon absorption
- `sigma_scat_alpha(α, ω, B)` — Thomson scattering (mode-dependent)

---

## What Needs Doing Next (Priority Order)

### 1. Fix Magnetic Scattering Separation (THE BLOCKER)
See "Immediate Blocker" above. This is ~30 lines of code in `magnetic_modes.jl` + `magnetic_atmosphere.jl`. After fixing, test:
```bash
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
# B=0 should match non-magnetic
r0 = solve_magnetic_atmosphere(1e6, 2e14, 0.0, 0.0, gaunt; K=50, M=8, N=100, tol=1e-4, verbose=true)
# B=10^12 should converge with F/σT⁴ ≈ 1
r1 = solve_magnetic_atmosphere(1e6, 2e14, 1e12, π/4, gaunt; K=50, M=8, N=100, tol=1e-4, verbose=true)'
```

### 2. Verify Magnetic Spectra Against Suleimanov+(2009)
Compare X-mode and O-mode emergent spectra at B=10¹² G, T_eff=10⁶ K against Suleimanov, Potekhin & Werner (2009) A&A 500, 891, Figs 4-6. The paper is NOT in `refs/` — download it.

### 3. Build SpectralImageCube v2
Replace the modified-blackbody placeholder in the rendering pipeline with real atmosphere spectra:
- Pre-compute atmosphere grid over (T_eff, log g, B, cos θ_e)
- Store as HDF5 with provenance metadata
- `SpectralImageCube` struct: see `NEUTRON_STAR_CLAUDE_MD_v2.md` for the full spec
- Re-render the tracer bullet with real atmosphere

### 4. (Optional) Frequency-Adaptive Depth Grid
Current ~5-10% spectrum error at peak is due to our fixed log(y) grid. McPHAC uses frequency-adaptive grids (`GetColumnsNu.c`):
- For each frequency k, restrict depth grid to τ < 80
- Log-spaced within that range (NDEPTHSNU points per frequency)
- Cubic spline interpolation of Eddington factors back to common grid
- Would improve accuracy from ~7% to <1%

### 5. (Phase 4) Extensions
In order: OCP Monte Carlo, partial ionisation, vacuum QED mode conversion, Kerr ray tracing, magnetosphere volume rendering. See `NEUTRON_STAR_CLAUDE_MD_v2.md` Phase 4.

---

## Architecture

### Module Dependency Order (NeutronStar.jl)
```
PhysicalConstants → BSkEOS → TOVSolver → DipoleModel
                  → GauntFactor → HydrogenOpacity → MagneticCoulomb → MagneticFF
                  → DielectricTensor → BlackbodyAtmosphere → MagneticModes
                  → AtmosphereStructure → FeautrierSolver → TemperatureCorrection
                  → RTAtmosphere → MagneticAtmosphere
                  → SchwarzschildTracer → CIE_sRGB → Renderer
```

### Key Files and What They Do

**Atmosphere (the active development area):**
| File | Purpose | Status |
|------|---------|--------|
| `src/atmosphere/feautrier.jl` | Block-tridiagonal Feautrier RT solver | Working, verified vs McPHAC |
| `src/atmosphere/temp_correction.jl` | Rybicki (1971) global temperature correction | Working, verified vs McPHAC |
| `src/atmosphere/rt_atmosphere.jl` | Non-magnetic atmosphere driver: iterates Feautrier + Rybicki | Working |
| `src/atmosphere/magnetic_atmosphere.jl` | Two-mode magnetic atmosphere driver | Framework done, B>0 blocked |
| `src/atmosphere/atm_structure.jl` | Column structure, Eddington T profile, opacity computation | Working |
| `src/atmosphere/blackbody.jl` | Planck function B_ν and modified blackbody | Working |

**Opacities:**
| File | Purpose | Status |
|------|---------|--------|
| `src/opacity/hydrogen_ff.jl` | B=0 free-free opacity (Haakonsen Eq. 12) + Thomson | Working |
| `src/opacity/gaunt_ff.jl` | Gaunt factor interpolation (Sutherland 1998 table) | Working |
| `src/opacity/magnetic_ff.jl` | Magnetic free-free σ_α cross-sections (P&C 2003) | Working |
| `src/opacity/magnetic_modes.jl` | Mode decomposition κ_j, angle-averaged K_∥/K_⊥ | Working, needs absorption/scattering split |
| `src/opacity/dielectric_tensor.jl` | Cold-plasma Stix parameters, polarization weights | Working |
| `src/opacity/coulomb_magnetic.jl` | Magnetic Coulomb logarithm (Eq. 44) | Working |

**Everything else (stable, don't touch unless needed):**
| File | Purpose |
|------|---------|
| `src/constants.jl` | CODATA physical constants |
| `src/eos/bsk_eos.jl` | BSk analytical EOS fits |
| `src/eos/tov.jl` | TOV integrator |
| `src/surface/dipole.jl` | Dipole B-field + Greenstein-Hartke T map |
| `src/geodesics/schwarzschild.jl` | GR ray tracing via elliptic integrals |
| `src/colorimetry/cie_srgb.jl` | CIE 1931 → XYZ → sRGB pipeline |
| `src/pipeline/render.jl` | End-to-end rendering orchestration |

### McPHAC Reference Data
```
refs/code/McPHAC/                  ← McPHAC source + binary
refs/code/McPHAC/OUT_T6_g2e14/    ← T_eff=10⁶, g=2×10¹⁴ (converged, N=400)
refs/code/McPHAC/OUT_T5.3_g1e14/  ← T_eff=10⁵·³, g=10¹⁴ (converged, N=400)
refs/code/McPHAC/OUT_T6.5_g2e14/  ← T_eff=10⁶·⁵, g=2×10¹⁴ (converged, N=400)
refs/code/McPHAC/gffgu.dat        ← Sutherland Gaunt factor table (used by both codes)
```
McPHAC output format: `Temperatures.{N}.{iter}.dat` has columns `#d logy oldT newT deltaT deltaT/T avgf |avgf/sigmaf|`. Use `newT` from the last iteration.

---

## Critical Physics Conventions

1. **Per-mode Planck source = B_ν/2**: Unpolarized thermal emission splits equally between X and O modes. The Feautrier source, bottom BC, and Rybicki K_k all use B_ν/2 per mode.

2. **Feautrier surface BC = pure radiation**: No local source (Q=0) and no scattering at the surface (i=1). McPHAC: CalcK.c `K[0]=0`, CalcU.c `U[0]=0`. Outgoing radiation comes through tridiagonal coupling to deeper layers.

3. **Rybicki surface = zero coupling**: `U_k[1]=0`, `K_k[1]=0`. The surface T_k uses `h_ν + f/Δ_surf` with `Δ_surf = 0.5×(k₁+k₂)×(y₂-y₁)` (note the 0.5 factor, different from interior Δ).

4. **EOS**: ρ = m_p P / (2 k_B T) for fully ionized hydrogen (μ=0.5). McPHAC uses OPAL EOS which gives the same at these conditions.

5. **B̄ in Rybicki is GLOBAL**: When summing over two modes, B̄_i = Σ_j Σ_k (B_ν/2)×κ_j×b_k / Σ_j Σ_k (dBdT/2)×κ_j×b_k. Do NOT compute per-mode — this was a bug that caused divergence.

6. **Flux diagnostic**: F = 2π ∫₀¹ I(μ) μ dμ (not 4π). The `_bolometric_flux` function uses the correct 2π factor.

---

## Quick Test Commands

```bash
# Compile check
julia --project=. -e 'using NeutronStar; println("OK")'

# Non-magnetic atmosphere (should converge in ~10 iter, F/σT⁴ ≈ 1.0)
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
result = solve_atmosphere(1e6, 2e14, gaunt; K=50, M=8, N=200, max_iter=30, tol=1e-6, verbose=true)'

# Magnetic atmosphere B=0 (should match non-magnetic within 0.1%)
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
r = solve_magnetic_atmosphere(1e6, 2e14, 0.0, 0.0, gaunt; K=50, M=8, N=100, tol=1e-4, verbose=true)'

# McPHAC reference run
cd refs/code/McPHAC && ./McPHAC 6.0 2e14 -9 2.5 100 4 100 4 8 50 1e-4 50 1 80 0.265
```

---

## Bugs Fixed This Session (for context, don't re-introduce)

| Bug | File | Root Cause | Fix |
|-----|------|-----------|-----|
| Surface Δ missing 0.5 | temp_correction.jl:206 | Paper convention differs at surface vs interior | `0.5 * (k₁+k₂) * dy` |
| Surface U_k nonzero | temp_correction.jl:216 | Was computing interior formula | Set `U_k[1] = 0` |
| Surface K_k nonzero | temp_correction.jl:219 | Was computing interior formula | Set `K_k[1] = 0` |
| Feautrier surface source | feautrier.jl:194-208 | Had `Q=(1-ρ)B` + scattering at i=1 | Skip source+scattering at i=1 |
| EOS μ=1 not μ=0.5 | atm_structure.jl:90 | `ρ=m_p P/(k_BT)` → 2× density | `ρ=m_p P/(2k_BT)` |
| Flux diagnostic 4π→2π | rt_atmosphere.jl:152 | Wrong solid angle factor | `2π` not `4π` |

## References

- `refs/haakonsen_2012_mcphac.pdf` — McPHAC paper (Appendix A = Rybicki method)
- `refs/potekhin_chabrier_2003_ff_opacity.pdf` — Magnetic opacities (Eqs. 25-52)
- `refs/potekhin_tables/hmagtab.txt` — Table format specification
- `NEUTRON_STAR_CLAUDE_MD_v2.md` — Master project specification (v2)
- **MISSING**: Suleimanov+ (2009) A&A 500:891 — need to download for magnetic verification
- **WARNING**: `refs/potekhin_chabrier_ho_2014_opacities.pdf` is the WRONG PAPER (Boev & Kovalev 2014 about exciton BEC). Needs replacement.
