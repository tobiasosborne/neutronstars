# Verification Log

## EOS (BSk21, BSk20, BSk19)

### Monotonicity
- **Test:** P(ρ) monotonically increasing for log ρ ∈ [6, 15.5]
- **Result:** PASS (all three EOS)
- **Source:** Potekhin et al. (2013) Fig. 2

### Maximum masses
| EOS | M_max computed | M_max published | Error | Status |
|-----|---------------|-----------------|-------|--------|
| BSk21 | 2.2735 M_☉ | 2.27 M_☉ | 0.15% | PASS |
| BSk20 | 2.1640 M_☉ | 2.16 M_☉ | 0.18% | PASS |
| BSk19 | 1.8598 M_☉ | 1.86 M_☉ | 0.01% | PASS |

**Source:** Potekhin et al. (2013) Section 6.

### Canonical 1.4 M_☉ NS (BSk21)
- M = 1.4000 M_☉, R = 12.59 km, ρ_c = 7.30×10¹⁴ g/cm³
- g_surface = 1.43×10¹⁴ cm/s², u = 0.3286
- **R error vs ~12.6 km from Fig. 8: 0.1%** → PASS

### Newtonian limit
- ρ_c = 5×10¹³: u = 0.0016 (≪ 1) → PASS

---

## Ray Tracer (Schwarzschild / Beloborodov 2002)

### Visible fraction
| u = r_g/R | f_vis computed | f_vis expected | Status |
|-----------|---------------|----------------|--------|
| 0 (flat)  | 0.5000        | 0.5000         | PASS   |
| 1/3 (R=3r_g) | 0.7500   | 0.7500         | PASS   |

**Source:** Beloborodov (2002) Eq. in Sect. 3.

### Flat-space limit
- M → 0: f_vis = 0.5001, u = 0.00025 → PASS

### Emission angles
- All cos α ∈ [0, 1] for hit rays → PASS

---

## Colorimetry

### Solar blackbody (5778K)
- Chromaticity (x, y) = (0.3264, 0.3357)
- Published Planckian at 5778K: ~(0.329, 0.339)
- Error: x = 0.8%, y = 1.0% → PASS (<2%)
- Normalized sRGB: (255, 252, 245) — warm white → PASS (R ≥ G ≥ B)

**Source:** CIE 015:2004 via CVRL data.

---

## Magnetic Opacity

### Vacuum-Polarized Normal Modes
- **Test:** Potekhin/Lai/Chabrier vacuum coefficient fits at `B=1e14 G`, mode parameters at `E=1 keV`, `ρ=1 g cm^-3`, `θ_B=45°`, and B-frame cyclic weights.
- **Source:** Potekhin, Lai & Chabrier (2004) vacuum coefficient fits and normal-mode equations; van Adelsberg & Lai (2006) B-frame cyclic component formulae.
- **Result:** PASS. Regression targets include `β=19.411715027`, `K1=-0.0257374555`, `K2=38.84916751`, and the expected `[|e_-|², |e_0|², |e_+|²]` weights for both modes.

### Suleimanov et al. (2009) Fig. 2
- **Test:** Digitized Fig. 2 and compared component opacities at `T=7e6 K`, `ρ=30 g cm^-3`, `B=1e14 G`, `θ_B={5,45,90} deg`.
- **Source:** Suleimanov, Potekhin & Werner (2009) Fig. 2; Potekhin & Chabrier (2003) Eq. 44a-e and Eq. 52.
- **Result:** PASS for the red-green free-free gate outside resonance-marker bands.
- `θ=5°`, magenta/free-free, mode 1: median abs residual `0.141 dex`.
- `θ=90°`, magenta/free-free, mode 1: median abs residual `0.081 dex`.
- **Artifacts:** `verification/figures/suleimanov_2009_fig2/opacity_comparison.png`, `verification/data/suleimanov_2009_fig2/opacity_comparison_metrics.json`.
- **Caveat:** This does not validate vacuum-resonance/mode-conversion physics; those bands are intentionally excluded from the metric.

### Suleimanov Fig. 2 Clean Branch View
- **Test:** Mask resonance-marker bands and tall annotation segments, split each color into upper/lower branches per x-column, and compare each branch against each computed mode without nearest-mode reassignment hiding the residuals.
- **Result:** FAIL for oblique-angle full-figure agreement.
- `θ=90°`: good agreement; free-free branch median abs residuals are about `0.08 dex`.
- `θ=5°`: mixed; best branch median abs residuals are about `0.4-0.5 dex`.
- `θ=45°`: major mismatch remains; upper magenta/blue branches show about `2.7-3.0 dex` median abs residual.
- **Artifacts:** `verification/figures/suleimanov_2009_fig2/clean_branch_validation.png`, `verification/data/suleimanov_2009_fig2/clean_branch_metrics.json`.
- **Next:** Isolate fixed-energy `θ_B=45°` intermediate quantities (`β`, `K_j`, `K_z,j`, cyclic weights, component cross-sections) against source equations.

---

## End-to-End Tracer Bullet (128×128)
- Parameters: M=1.4 M_☉, R=12 km, B=10¹² G, T_pole=10⁶ K, T_eq=5×10⁵ K
- True colour: uniform blue-white (physically correct for 10⁶ K BB)
- False colour: shows hotspot structure, full dynamic range
- 71.3% of pixels hit the star (consistent with f_vis = 0.763 for u = 0.345)
- Per-pixel spectral data saved

**Status:** Phase 2 tracer bullet OPERATIONAL.
