# Session Handoff — 2026-03-24

## Project Overview

Building a physically exact neutron star visualisation pipeline in Julia. Every equation must be traceable to a locally-stored peer-reviewed PDF. The pipeline goes: TOV(M,R,g) → surface model T(θ,B) → radiative transfer → GR ray tracing → colorimetry → rendered image.

## What Exists and Works

### Phase 2: Tracer Bullet — COMPLETE, ALL VERIFIED
- **BSk EOS** (`src/eos/bsk_eos.jl`): BSk19/20/21 from Potekhin+ (2013) Eq. 3. M_max within 0.2% of published values.
- **TOV solver** (`src/eos/tov.jl`): RK4 adaptive. Canonical 1.4 M_☉ → R=12.59 km.
- **Schwarzschild ray tracer** (`src/geodesics/schwarzschild.jl`): Beloborodov (2002) cosine relation. Visible fraction exact. Flat-space limit correct.
- **Surface model** (`src/surface/dipole.jl`): Greenstein-Hartke T(θ_B), dipole B(θ_B).
- **Colorimetry** (`src/colorimetry/cie_srgb.jl`): CIE 1931 → XYZ → sRGB → Reinhard tonemapping. Solar BB chromaticity within 1%.
- **Pipeline** (`src/pipeline/render.jl`): End-to-end 512×512 renders in `output/`.

### Phase 3.1: B=0 RT Solver — FUNCTIONAL, NEEDS TEMP CORRECTION
- **Gaunt factors** (`src/opacity/gaunt_ff.jl`): Sutherland (1998) table, 81×41 grid, bilinear interpolation. Verified.
- **B=0 opacities** (`src/opacity/hydrogen_ff.jl`): Free-free κ_ff (Haakonsen Eq. 12) + Thomson σ. Matches Kramers ±30%.
- **Atmosphere structure** (`src/atmosphere/atm_structure.jl`): Hydrostatic P=g_s y, Eddington T(y), ideal gas EOS. Auto-extends y_max until τ≥80.
- **Feautrier solver** (`src/atmosphere/feautrier.jl`): Block-tridiagonal with Auer (1976) discretisation. J→B_ν at depth (exact). f→1/3 at depth (exact). All emergent I > 0. Anisotropic Thomson scattering (Haakonsen Eq. 16). **Fixed bug:** isotropic scattering factor was 2× too large (cross=1.0 should have been cross=0.5); corrected.
- **RT driver** (`src/atmosphere/rt_atmosphere.jl`): Single-pass (no temperature iteration). Spectral hardening I/B_ν ≈ 1.5 near peak. Flux ratio F/σT⁴ ≈ 2.8.
- **Known limitation:** Temperature correction (`src/atmosphere/temp_correction.jl`) uses Unsöld-Lucy local correction which DIVERGES for scattering-dominated atmospheres. The Rybicki method (Haakonsen Appendix A, Eqs. A6-A33) is needed. This is a global tridiagonal system coupling all depths and frequencies. Not yet implemented.

### Phase 3.2: Magnetic Opacities — IN PROGRESS
- **Magnetic cross-sections** (`src/opacity/magnetic_ff.jl`): Three polarisation modes (α=-1,0,+1). Eq. 37-52 from Potekhin & Chabrier (2003). Electron cyclotron resonance, proton cyclotron, combined damping.
- **Magnetic Coulomb logarithm** (`src/opacity/coulomb_magnetic.jl`): Full Landau-level sum (Eq. 44a-e). At β_e=134 (B=10¹², T=10⁶), gives Λ_mag/Λ_classical = 0.053 (19× reduction). Uses QuadGK adaptive integration + SpecialFunctions besselk.
- **Mode decomposition** (`src/opacity/magnetic_modes.jl`): Polarisation vectors (simplified cold-plasma), angle-averaged κ_∥/κ_⊥ (Eq. 29-30), Rosseland frequency mean.
- **Potekhin tables downloaded** to `refs/potekhin_tables/` (47 .dat files, B=10^{10.5}-10^{15} G).

## Current Accuracy vs Potekhin Tables

At B=10¹² G, T=10⁶ K (fully ionised regime):

```
lg(ρ)  lg(K∥_ours)  lg(K∥_table)  Δ∥     lg(K⊥_ours)  lg(K⊥_table)  Δ⊥
-5.0    -1.37        -2.71        +1.34   -1.03         -3.03        +2.00
-3.0    -1.35        -2.29        +0.95   -1.01         -2.56        +1.55
-1.0    -1.15        -0.81        -0.35   -0.86         -1.08        +0.22
 0.0    -0.79        +0.06        -0.85   -0.56         -0.21        -0.35
 1.0    -0.10        +0.89        -0.98   +0.09         +0.65        -0.56
```

**Pattern:** Too high at low ρ (polarisation weights wrong), too low at high ρ (overcorrected). The K_⊥ is consistently closer than K_∥, suggesting the angle averaging weights cos²θ and sin³θ differently from correct.

## What Needs Doing Next (Priority Order)

### 1. Full Dielectric Tensor for Polarisation Vectors (~biggest remaining fix)

**File:** Create `src/opacity/dielectric_tensor.jl` or rewrite `polarisation_weights()` in `magnetic_modes.jl`.

**Problem:** Current `polarisation_weights()` uses a crude cold-plasma approximation that doesn't properly split opacity between the extraordinary (j=1) and ordinary (j=2) modes. At low density where the plasma frequency is small, the mode decomposition is especially wrong.

**Solution:** Implement the full cold-plasma dielectric tensor (Ginzburg 1970 §10):
```
ε_xx = 1 - ω_pe²ω/[ω(ω²-ω_ce²)] - ω_pi²ω/[ω(ω²-ω_cp²)]
ε_xy = i[ω_pe²ω_ce/(ω(ω²-ω_ce²)) - ω_pi²ω_cp/(ω(ω²-ω_cp²))]
ε_zz = 1 - ω_pe²/ω² - ω_pi²/ω²
```

Then compute the mode parameter q (Eq. 25) and derive |e_{j,α}|² properly. The paper references Shafranov (1967) and Ho & Lai (2003) for the polarisation vector formulae.

**Key equation:** Eq. 25 of P&C 2003:
```
q + ip = [ε_yy - ε_xx cos²θ_B + ε_zz sin²θ_B] / [2i(ε_xy cosθ_B + ε_yz sinθ_B)]
```

From q and p, the polarisation weights follow from the normal mode eigenvectors of the wave equation in the magnetised plasma.

### 2. Non-Born Gaunt Factor Correction (~0.2 dex)

**Where:** In `magnetic_ff.jl`, function `nu_ff_alpha()`.

**What:** The Born approximation Coulomb logarithm (Eq. 44) underestimates the true free-free absorption at low frequencies. Correct by multiplying Λ_α^ff by the ratio g̃^ff/g̃_Born, where g̃^ff is the non-Born Gaunt factor from Hummer (1988) Padé formula. The paper says (top of p. 966): "we multiply Λ_α^ff by the ratio g̃^ff/g̃^ff_Born where g̃^ff_Born = (√3/π)Λ_cl^ff".

### 3. Better Numerical Resolution

**Frequency grid:** Current 60 log-spaced points may miss cyclotron resonance structure. Need ≥200 points, with additional points clustered near ω_ce and ω_cp.

**Angle grid:** Current 20-point midpoint rule in `kappa_parallel_mono`/`kappa_perp_mono`. Need ≥40 Gauss-Legendre points, with clustering near mode crossings.

### 4. Rybicki Temperature Correction (for Phase 3.1 McPHAC match)

**File:** `src/atmosphere/temp_correction.jl` needs complete rewrite.

**What:** Implement Haakonsen Appendix A (Eqs. A6-A33). This is a Rybicki (1971) method: for each frequency k, build a tridiagonal matrix T_k and coupling vectors U_k, V_k. Accumulate dense N×N matrix W = -I + Σ_k V_k T_k⁻¹ U_k. Solve W × ΔT = Q. Back-substitute per frequency.

**Why current method fails:** The Unsöld-Lucy local correction ΔT = -∫(J-B)κdν / ∫(dB/dT)κdν is unstable in scattering-dominated atmospheres because J is determined non-locally. The correction keeps pushing T in one direction without convergence. The Rybicki method accounts for non-local radiative coupling.

### 5. Partial Ionisation (for table match at T < 10⁶ K)

**Cannot be done from the paper alone.** The bound-bound and bound-free cross-sections come from precalculated quantum-mechanical tables (Forster+ 1984, Potekhin & Pavlov 1997) that are not published. Options:
- Email Potekhin (palex@astro.ioffe.ru) for the Fortran code/data
- Use the tables themselves for the partially ionised regime (verification only, per project rules)
- Implement the ionisation equilibrium from the EOS equations (Sect. 2.2, Eq. 19-23) and use the fully-ionised opacity model, accepting the error at low T

## File Map

```
src/
├── NeutronStar.jl              ← Main module. LOAD ORDER MATTERS (see below)
├── constants.jl                ← CODATA 2018, CGS units
├── eos/
│   ├── bsk_eos.jl              ← BSk19/20/21 P(ρ) from Potekhin+ 2013
│   └── tov.jl                  ← TOV solver, RK4 adaptive
├── surface/
│   └── dipole.jl               ← Greenstein-Hartke T(θ), dipole B(θ)
├── opacity/
│   ├── gaunt_ff.jl             ← Sutherland Gaunt factor table
│   ├── hydrogen_ff.jl          ← B=0 free-free + Thomson
│   ├── coulomb_magnetic.jl     ← Magnetic Coulomb log Λ_α^ff (Eq. 44)
│   ├── magnetic_ff.jl          ← Magnetic σ_α^ff/pp/scat (Eq. 37-52)
│   └── magnetic_modes.jl       ← Mode decomposition + Rosseland K_∥/K_⊥
├── atmosphere/
│   ├── blackbody.jl            ← Planck function, modified BB (Phase 2)
│   ├── atm_structure.jl        ← Column struct, hydrostatic eq, Eddington T(y)
│   ├── feautrier.jl            ← Block-tridiagonal Feautrier RT solver
│   ├── temp_correction.jl      ← Unsöld-Lucy (BROKEN for scattering; needs Rybicki)
│   └── rt_atmosphere.jl        ← Single-pass RT driver
├── geodesics/
│   └── schwarzschild.jl        ← Beloborodov (2002) ray tracer
├── colorimetry/
│   └── cie_srgb.jl             ← CIE 1931 → XYZ → sRGB → Reinhard
└── pipeline/
    └── render.jl               ← End-to-end rendering pipeline
```

**Module load order in NeutronStar.jl** (dependencies flow downward):
1. PhysicalConstants
2. BSkEOS → TOVSolver
3. DipoleModel
4. GauntFactor → HydrogenOpacity → MagneticCoulomb → MagneticFF
5. BlackbodyAtmosphere → MagneticModes
6. AtmosphereStructure → FeautrierSolver → TemperatureCorrection → RTAtmosphere
7. SchwarzschildTracer → CIE_sRGB
8. Renderer

**Critical:** `MagneticModes` depends on `BlackbodyAtmosphere` (uses `planck_Bnu` and `dBnu_dT`). `MagneticFF` depends on `MagneticCoulomb`. If you add new modules, check the dependency chain.

## References

All papers in `refs/`:
- **EOS/TOV:** `potekhin_2013_bsk_eos.pdf` — Eq. 3, Table 2
- **Magnetic opacity:** `potekhin_chabrier_2003_ff_opacity.pdf` — Eqs. 27-52, Appendix A+B (THIS IS THE KEY PAPER for the table reproduction work)
- **B=0 atmosphere:** `haakonsen_2012_mcphac.pdf` — Feautrier method, Eqs. 2-22, Appendix A
- **Ray tracing:** `beloborodov_2002_ray_tracing.pdf` — cosine relation Eq. 1
- **Surface temp:** `greenstein_hartke_1983_surface_T.pdf`
- **Potekhin tables:** `refs/potekhin_tables/hmag*.dat` — verification targets (47 files)
- **McPHAC code:** `refs/code/McPHAC/` — C source, structural reference for B=0 solver
- **af proof tree:** `proof/` — 38 nodes, 19 external references, all linked to local PDFs

## Quick Test Commands

```bash
# Compile check
cd /home/tobias/Projects/neutronstars
julia --project=. -e 'using NeutronStar; println("OK")'

# EOS verification (should all PASS)
julia --project=. verification/eos_checks.jl

# Magnetic opacity vs Potekhin table
julia --project=. -e '
using NeutronStar; using Printf; using NeutronStar.PhysicalConstants: k_B, h
B=1e12; T=1e6
ν_grid=[10.0^logν for logν in range(log10(0.05*k_B*T/h), log10(120*k_B*T/h), length=60)]
for (lgR,lgK0,lgK1) in [(-3.0,-2.293,-2.563),(0.0,0.061,-0.206),(1.0,0.887,0.652)]
    ρ=10.0^lgR; Kp,Kt=rosseland_magnetic(B,T,ρ,ν_grid)
    @printf("lgρ=%+.0f: K∥=%.2f(tab %.2f) K⊥=%.2f(tab %.2f)\n", lgR, log10(max(Kp,1e-30)), lgK0, log10(max(Kt,1e-30)), lgK1)
end'

# Render tracer bullet (takes ~2 min at 512×512)
julia --project=. -e '
using NeutronStar; using NeutronStar.Renderer: save_ppm
p=NSParams(1.4,12.0,1e12,1e6,5e5,π/6,π/3,100.0,1.5)
r=render_neutron_star(p,128;n_ν=80); save_ppm(r,"output/test")'
```
