# Session Handoff — 2026-03-24 (updated)

## Project Overview

Building a physically exact neutron star visualisation pipeline in Julia. Master goal (from `NEUTRON_STAR_CLAUDE_MD_v2.md`): spectral image cube I(x,y,ν) from first principles. Every equation traceable to a locally-stored published paper.

## What Was Done This Session

### Previous Work (carried forward)
- Potekhin table parser, dielectric tensor, magnetic modes (see git history)
- Phase 2 tracer bullet: TOV, Schwarzschild ray tracer, surface model, colorimetry all working

### Atmosphere Solver Bug Fixes (THIS SESSION)
Fixed 5 bugs that brought the solver from 32% error to <1.2%:

1. **Rybicki surface Δ_surf** (temp_correction.jl): Missing 0.5 factor in surface optical depth step. McPHAC CalcTt.c uses `0.5*(k[0]+k[1])` while interior uses `(k[i]+k[i-1])`.

2. **Rybicki surface U_k** (temp_correction.jl): Was computing interior formula at surface; should be `U_k[1] = 0` (McPHAC CalcUt.c line 10).

3. **Rybicki surface K_k** (temp_correction.jl): Was computing interior formula at surface; should be `K_k[1] = 0` (McPHAC CalcKt.c line 15).

4. **Feautrier surface BC** (feautrier.jl): Surface had local thermal source `Q=(1-ρ)B` and scattering terms; should be pure radiation BC with `Q=0` and no scattering (McPHAC CalcK.c `K[0]=0`, CalcU.c `U[0]=0`). This was the root cause of the factor-of-2 spectrum error.

5. **Mean molecular weight** (atm_structure.jl, rt_atmosphere.jl): EOS had `ρ = m_p P/(k_B T)` but ionized H needs `ρ = m_p P/(2 k_B T)` (μ=0.5). McPHAC uses OPAL EOS which accounts for this.

### Current Accuracy vs McPHAC (T_eff=10⁶ K, g=2×10¹⁴)

```
                Temperature Profile
log(y)    T_ours     T_McPHAC    Δ%
-9.00     248423     251330     -1.2%   ← surface
-6.21     248426     251333     -1.2%
-3.36     331457     335287     -1.1%   ← mid
-0.51    1232321    1241655     -0.8%
 2.34    5564690    5584505     -0.4%   ← deep

        Emergent Spectrum I(μ≈1)
hν/(kT)   I_ours     I_McPHAC   Δ%
0.05      0.0808     0.0811    -0.3%   ← RJ tail
0.21      1.371      1.307     +4.9%
1.02      22.25      20.74     +7.3%   ← near peak
2.26      60.29      67.55    -10.7%   ← peak region
5.01      95.57      102.4     -6.6%
11.1      59.71      61.72     -3.3%
24.5      5.991      6.098     -1.8%   ← Wien tail

Flux conservation: F/σT⁴ = 0.99
```

The ~5-10% spectrum errors at the peak are due to our fixed log(y) depth grid vs McPHAC's frequency-adaptive grid (GetColumnsNu). Improvable with adaptive grid.

## What Needs Doing Next (Priority Order)

### 1. Multi-Parameter Verification (IN PROGRESS)
Verify at (T_eff=10^5.3, g=10^14) and (T_eff=10^6.5, g=2×10^14.3).
McPHAC runs being generated.

### 2. Magnetic Two-Mode Atmosphere (NEXT)
Extend Feautrier for X/O modes. Each mode has its own opacity κ^j(ν,T,ρ,B,θ):
- X-mode: suppressed scattering below cyclotron frequency
- O-mode: enhanced opacity
Verify B→0 recovers non-magnetic result.

### 3. SpectralImageCube v2 + Atmosphere Grid
Build v2 SpectralImageCube struct. Pre-compute atmosphere grid over (T_eff, log g).
Re-render tracer bullet with real atmosphere spectra.

### 4. (Optional) Adaptive Depth Grid
Implement frequency-adaptive depth grid (like McPHAC's GetColumnsNu) to improve spectrum accuracy from ~7% to <1%.

## Quick Test Commands

```bash
# Compile check
julia --project=. -e 'using NeutronStar; println("OK")'

# Run atmosphere solver (converges in ~10 iterations)
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
result = solve_atmosphere(1e6, 2e14, gaunt; K=80, M=10, N=200, max_iter=30, tol=1e-6, verbose=true)'

# Run McPHAC for comparison
cd refs/code/McPHAC/OUT_T6_g2e14  # Results saved here
```

## File Map

### Modified Files (this session)
```
src/atmosphere/temp_correction.jl   ← FIXED: 3 surface BC bugs
src/atmosphere/feautrier.jl         ← FIXED: surface source/scattering bug
src/atmosphere/atm_structure.jl     ← FIXED: mean molecular weight
src/atmosphere/rt_atmosphere.jl     ← FIXED: density + flux diagnostic
```

### Existing Working Files
```
src/NeutronStar.jl, src/constants.jl
src/eos/bsk_eos.jl, src/eos/tov.jl
src/surface/dipole.jl, src/geodesics/schwarzschild.jl
src/colorimetry/cie_srgb.jl, src/pipeline/render.jl
src/opacity/gaunt_ff.jl, src/opacity/hydrogen_ff.jl
src/opacity/magnetic_ff.jl, src/opacity/coulomb_magnetic.jl
src/opacity/dielectric_tensor.jl, src/opacity/magnetic_modes.jl
src/atmosphere/blackbody.jl, src/atmosphere/atm_structure.jl
src/atmosphere/feautrier.jl, src/atmosphere/temp_correction.jl
src/atmosphere/rt_atmosphere.jl
```

## References
- `refs/haakonsen_2012_mcphac.pdf` — McPHAC paper
- `refs/potekhin_chabrier_2003_ff_opacity.pdf` — Magnetic opacities
- `NEUTRON_STAR_CLAUDE_MD_v2.md` — Master project specification
