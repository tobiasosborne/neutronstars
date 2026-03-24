# Session Handoff — 2026-03-24

## Project Overview

Building a physically exact neutron star visualisation pipeline in Julia. Master goal (from `NEUTRON_STAR_CLAUDE_MD_v2.md`): spectral image cube I(x,y,ν) from first principles. Every equation traceable to a locally-stored published paper.

## What Was Done This Session

### 1. Potekhin Table Parser (DONE)
- `src/eos/potekhin_table_reader.jl`: Parses all 47 .dat files (66,390 data points)
- Handles overlapping B-value files via `merge_tables()`
- `verification/potekhin_table_comparison.jl`: Automated comparison harness

### 2. Dielectric Tensor for Normal Mode Polarization (DONE)
- `src/opacity/dielectric_tensor.jl`: Full cold-plasma Stix parameters (S, D, P), q-based polarization vectors (avoids catastrophic cancellation at low density)
- `src/opacity/magnetic_modes.jl`: REWRITTEN — fixed two bugs:
  1. Polarization weights from full wave equation (replaces ad-hoc Eq. 26 approximation)
  2. Angle-average each mode separately THEN combine (P&C 2003 Eq. 30)
- Verified: fully-ionized low-density opacity error reduced from 1.34 dex to 0.17 dex

### 3. Rybicki Temperature Correction (IMPLEMENTED, NEEDS DEBUGGING)
- `src/atmosphere/temp_correction.jl`: REWRITTEN with Rybicki (1971) method from Haakonsen (2012) Appendix A
- **Algorithm**: For each frequency k, build tridiagonal T_k, solve N column-wise for T_k⁻¹×diag(U_k), accumulate dense W matrix, solve (W+I)J̄ = rhs
- **Converges** in 9 iterations (max|ΔT/T| < 10⁻⁴)
- **But**: temperature profile doesn't match McPHAC — 15% error at spectral peak

### 4. RT Atmosphere Driver (REWRITTEN)
- `src/atmosphere/rt_atmosphere.jl`: Full iteration loop with damped Rybicki correction
- Builds initial Eddington T(y), iterates RT + temperature correction until convergence
- McPHAC compiled and running as verification target (`refs/code/McPHAC/OUT/`)

### 5. Project Reframing
- Read and analyzed `NEUTRON_STAR_CLAUDE_MD_v2.md`
- **Key insight**: Potekhin tables are VERIFICATION targets, not production targets
- **Key decision**: Fully ionized H only for Phase 3. Partial ionization = Phase 4.
- Plan updated in `.claude/plans/giggly-tickling-quasar.md`

## Current Accuracy vs McPHAC (T_eff=10⁶ K, g=2×10¹⁴)

```
                    Temperature Profile
log(y)    T_ours     T_McPHAC    ratio
-9.00     171605     251329      0.683  ← surface too cold
-4.54     170778     256606      0.666
-2.25     458752     529384      0.867
 0.03    1751479    1650144      1.061  ← deep 6% too hot
 2.32    5930456    5584505      1.062

              Spectrum (F_ν*ν/σT⁴)
hν/(kT)   McPHAC       Ours         Δ%
0.13      7.85e-05    8.10e-05     +3.2%   ← good
0.87      1.60e-02    1.31e-02    -18.4%   ← peak region
2.26      1.35e-01    1.17e-01    -13.3%
5.87      5.00e-01    5.52e-01    +10.4%   ← Wien tail excess
15.2      3.89e-01    6.37e-01    +64.1%   ← deep T too hot
```

## Known Bug: Rybicki Temperature Profile Shape

The Rybicki correction converges but to the wrong temperature profile. Surface is 32% too cold, deep atmosphere 6% too hot. This distorts the spectrum.

**Root cause hypothesis**: The discretization of the temperature correction equation (Haakonsen A11) may have remaining issues in how the surface boundary condition (A28-A30) couples to the interior Rybicki system. The Δ convention (A12-A14 — paper omits 1/2 factor) was already fixed.

**Debug leads**:
1. Compare the T_k matrix elements against McPHAC's C code (`DeltaT.c`, `CalcU.c`, `CalcV.c`)
2. Check if h_ν (flux Eddington factor at surface) used in A30 matches McPHAC
3. Try the simpler approach: use Unsöld-Lucy for the first few iterations (to get close), then switch to Rybicki (the U-L divergence only matters at convergence, not early iterations)
4. Check the normalization of dB/dT in the Rybicki U_k and V_k vectors

## File Map

### New/Modified Files (this session)
```
src/eos/potekhin_table_reader.jl         ← NEW: Potekhin table parser
src/opacity/dielectric_tensor.jl         ← NEW: Full cold-plasma normal modes
src/opacity/magnetic_modes.jl            ← REWRITTEN: Proper mode decomposition
src/atmosphere/temp_correction.jl        ← REWRITTEN: Rybicki method
src/atmosphere/rt_atmosphere.jl          ← REWRITTEN: Iterative RT driver
src/NeutronStar.jl                       ← MODIFIED: Added DielectricTensor module
verification/potekhin_table_comparison.jl ← NEW: Comparison harness
```

### Existing Working Files (Phase 2, unchanged)
```
src/constants.jl, src/eos/bsk_eos.jl, src/eos/tov.jl
src/surface/dipole.jl, src/geodesics/schwarzschild.jl
src/colorimetry/cie_srgb.jl, src/pipeline/render.jl
src/opacity/gaunt_ff.jl, src/opacity/hydrogen_ff.jl
src/opacity/magnetic_ff.jl, src/opacity/coulomb_magnetic.jl
src/atmosphere/blackbody.jl, src/atmosphere/atm_structure.jl
src/atmosphere/feautrier.jl
```

## What Needs Doing Next (Priority Order)

### 1. Fix Rybicki Temperature Profile (~biggest remaining fix)
Debug the temperature correction to match McPHAC within 1%. See "Debug leads" above. The most productive approach is probably to compare the actual T_k matrix elements line-by-line against McPHAC's C code.

### 2. Verify B=0 Atmosphere Against McPHAC
Once T profile matches, verify emergent spectrum at 3+ (T_eff, log g) combinations: (10⁵·³, 14.0), (10⁶, 14.3), (10⁶·⁵, 14.3). Target: <1%.

### 3. Magnetic Two-Mode Atmosphere
Extend Feautrier for two decoupled polarization modes. Use existing dielectric tensor + magnetic opacity modules. Verify B→0 recovers non-magnetic result.

### 4. SpectralImageCube + Atmosphere Grid
Build the v2 SpectralImageCube struct. Pre-compute atmosphere grid. Re-render tracer bullet.

## Quick Test Commands

```bash
# Compile check
julia --project=. -e 'using NeutronStar; println("OK")'

# Run atmosphere solver
julia --project=. -e '
using NeutronStar; using NeutronStar.GauntFactor: load_gaunt_table
gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
result = solve_atmosphere(1e6, 2e14, gaunt; K=50, M=8, N=100, max_iter=30, tol=1e-4, verbose=true)'

# Run McPHAC for comparison
cd refs/code/McPHAC
./McPHAC 6.0 2e14 -9 2.5 100 4 100 4 8 50 1e-4 50 1 80 0.265
# Output in OUT/SurfaceFluxesNorm.400.6.dat, OUT/Fluxes.400.6.dat

# Magnetic opacity check
julia --project=. -e '
using NeutronStar; using NeutronStar.PhysicalConstants: k_B, h
B=1e12; T=1e6
ν_grid=[10.0^logν for logν in range(log10(0.05*k_B*T/h), log10(120*k_B*T/h), length=60)]
ρ=1e-3; Kp,Kt=rosseland_magnetic(B,T,ρ,ν_grid)
println("K∥=$(log10(Kp)) K⊥=$(log10(Kt))  [table: -2.29, -2.56]")'
```

## References

- `refs/haakonsen_2012_mcphac.pdf` — McPHAC paper, Appendix A = Rybicki temperature correction
- `refs/potekhin_chabrier_2003_ff_opacity.pdf` — Magnetic opacities (Eqs. 25-52)
- `refs/potekhin_tables/hmagtab.txt` — Table format specification
- `NEUTRON_STAR_CLAUDE_MD_v2.md` — Master project specification (v2)
- `.claude/plans/giggly-tickling-quasar.md` — Current implementation plan

## IMPORTANT NOTE
`refs/potekhin_chabrier_ho_2014_opacities.pdf` is the WRONG PAPER — it's Boev & Kovalev (2014) about exciton BEC. Needs to be replaced with the actual Potekhin, Chabrier & Ho (2014) A&A 572, A69.
