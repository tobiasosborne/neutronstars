# Session Handoff — 2026-03-24

## What Was Done

### Phase 1: Literature & Semantic Graph
- **19 papers downloaded** to `refs/` (all arXiv + Springer Ruder book + Reinhard tonemapping + Greenstein/Pechenick from ADS)
- **CIE 1931 colour matching functions** downloaded from CVRL
- **McPHAC source code** cloned to `refs/code/McPHAC/`
- **Potekhin opacity tables** downloaded to `refs/potekhin_tables/` (47 .dat files, B=10^{10.5}-10^{15} G)
- **af proof tree**: 38 nodes across 6 pipeline stages, 19 external references
- **Fixed 2 paper mismatches**: haakonsen_2012 (was math paper, now correct McPHAC), potekhin_chabrier_2003 (was galaxy cluster paper, now correct magnetic opacity paper)
- **4 books still need manual browser download** (Wiley Cloudflare blocks headless): Rybicki & Lightman, Shapiro & Teukolsky, Mészáros, Mihalas

### Phase 2: Tracer Bullet — COMPLETE
All verification passing:
- **BSk EOS** (3 models): M_max within 0.2% of Potekhin+ (2013)
- **TOV solver**: RK4 adaptive, canonical 1.4 M_☉ → R=12.59 km
- **Schwarzschild ray tracer**: Beloborodov (2002), visible fraction exact, flat-space limit correct
- **Surface model**: Greenstein-Hartke T(θ), dipole B(θ)
- **Colorimetry**: CIE → sRGB → Reinhard, solar BB chromaticity within 1%
- **512×512 render** saved in `output/`

### Phase 3: RT Solver — IN PROGRESS
7 new files implemented:
- `src/opacity/gaunt_ff.jl` — Sutherland Gaunt factor table (81×41 grid) ✓
- `src/opacity/hydrogen_ff.jl` — B=0 free-free + Thomson ✓ (matches Kramers ±30%)
- `src/opacity/magnetic_ff.jl` — Magnetic cross-sections (Eq. 37-52) ✓
- `src/opacity/magnetic_modes.jl` — Mode decomposition + Rosseland K_∥/K_⊥ ✓
- `src/atmosphere/atm_structure.jl` — Column structure + Eddington T(y) ✓
- `src/atmosphere/feautrier.jl` — Block-tridiagonal Feautrier solver ✓ (J→B at depth, f→1/3)
- `src/atmosphere/rt_atmosphere.jl` — Single-pass driver ✓ (spectral hardening I/B~1.5)

### Magnetic Opacity Status
First comparison against Potekhin tables (B=10¹², T=10⁶):
- **Right order of magnitude** — magnetic suppression working (~100× below B=0)
- **Still ~0.5-2 dex high** at low density — using classical Coulomb log instead of full magnetic (Eq. 44)
- K_∥/K_⊥ ratio has correct sign but too small

## What Needs Doing Next

### Immediate: Improve Magnetic Opacity (0.5-2 dex discrepancy)
1. **Implement magnetic Coulomb logarithm** (Eq. 44a-e) — replaces classical Λ_cl with Λ_α^ff that accounts for Landau quantization. This is the main source of the remaining discrepancy. Requires summing over Landau levels n and integrating modified Bessel functions.
2. **Fix polarisation vector decomposition** in `magnetic_modes.jl` — current simplified cold-plasma weights are approximate. Need full dielectric tensor treatment (Eq. 25).
3. **Verify at multiple B values** (10^{11}, 10^{12}, 10^{13}) and check B→0 recovery.

### Next: Iterative Temperature Correction
- The Unsöld-Lucy local correction (Eq. A5) diverges for scattering-dominated atmospheres
- Need full **Rybicki method** (Haakonsen Appendix A, Eqs. A6-A33): global tridiagonal system coupling all depths and frequencies
- This is required to get <1% McPHAC agreement

### Then: Pipeline Integration
- Add `atmosphere_model` field to `NSParams` (`:blackbody` or `:rt`)
- Modify `render.jl` to use RT atmosphere with angle-dependent emergent spectra
- Pre-compute atmosphere grid for efficiency

## File Map (src/)
```
constants.jl              — CODATA 2018 physical constants (CGS)
eos/bsk_eos.jl            — BSk19/20/21 EOS (Potekhin+ 2013 Eq. 3)
eos/tov.jl                — TOV solver (RK4 adaptive)
surface/dipole.jl         — Greenstein-Hartke T(θ) + dipole B(θ)
opacity/gaunt_ff.jl       — Sutherland Gaunt factor table
opacity/hydrogen_ff.jl    — B=0 free-free + Thomson
opacity/magnetic_ff.jl    — Magnetic free-free (P&C 2003 Eqs. 37-52)
opacity/magnetic_modes.jl — Mode decomposition + Rosseland K_∥/K_⊥
atmosphere/blackbody.jl   — Modified BB placeholder (Phase 2)
atmosphere/atm_structure.jl — Column structure + Eddington T(y)
atmosphere/feautrier.jl   — Feautrier block-tridiagonal RT solver
atmosphere/temp_correction.jl — Unsöld-Lucy correction (needs Rybicki upgrade)
atmosphere/rt_atmosphere.jl — RT driver + emergent_spectrum interface
geodesics/schwarzschild.jl — Beloborodov (2002) ray tracer
colorimetry/cie_srgb.jl   — CIE 1931 → XYZ → sRGB → Reinhard
pipeline/render.jl         — End-to-end rendering pipeline
NeutronStar.jl            — Main module (load order matters!)
```
