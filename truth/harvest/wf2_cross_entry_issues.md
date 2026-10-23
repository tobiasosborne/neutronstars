## harvest:colour_srgb skipped
- Tristimulus integrals X = k * integral(S(lambda) xbar(lambda) dlambda): not stated in the README or the sRGB text.
- CIE D65 spectral power distribution and CIE standard illuminant tables: not in local sources.
- XYZ normalisation (Y=1 for white, k constant): not stated by Stokes 1996.
- Rounding rule from 255.0*R' to an integer code value: not stated.
- IEC 61966-2-1 thresholds 0.0031308 / 0.04045 and the 7-digit matrix: not in refs/.
- D50 ambient white, digital-TV counts 16/235, BT.709 transfer-function details: kept out of scope.

## harvest:condensed skipped
- Potekhin (1999) eqs. 28/39 for the collision times tau_par, tau_perp: formulas not in refs/
- Medin & Lai (2006, 2007) cohesive-energy tables and vapour-density method: only the cited T_crit fit is in Potekhin 2012
- H-surface T_crit values (8e4, 5e5, 1e6 K at 1e13, 1e14, 5e14 G) and Fe upper limit Q_s <~ Z^(9/5) B12^(2/5) eV, T_crit <~ 10^5.5 B12^(2/5) K: kept only in the notes of rt.cond-tcrit-fit-fe, no separate entries (file size cap)
- Coulomb-coupling melting criterion (van Adelsberg 2005 sect. 2.1: Gamma = 0.227 Z^2 (rho1/A)^(1/3)/T6 < Gamma_m = 175; solid for T < 7e4 eta^(1/3) B14^(2/5) K for H, 4e6 eta^(1/3) B12^(2/5) K for Fe): not entered
- van Adelsberg 2005 eq. 13 (dielectric tensor rotated to the surface frame), eqs. 16, 18, 20 (normal-incidence R and index approximations near E_Bi, E_pe), and Appendix A/B: not entered as separate entries; eq. 16 is noted under rt.cond-free-ion-analytic-reflectivity
- Potekhin 2012 appendix A quartic mode polynomial (A.6-A.21), root selection Im(n)<=0 (A.23) and Fresnel boundary conditions (A.27-A.29), plus B.8-B.12 scalar products: only referenced in the reimplement field and notes
- Potekhin 2012 temperature-correction scheme for the last atmosphere point (Unsoeld-Lucy, eqs. 31-34): referenced in notes of rt.cond-inner-bc, not entered
- Collisional damping numerical values (Fig. 1 of van Adelsberg 2005) and all figure data: not in extracted text
- Kerr/Doppler or flux integration beyond F_nu = int J B_nu cos theta dOmega: only in notes of rt.cond-kirchhoff-law
- Universal constants hbar, e, c, m_e, k_B, m_u omitted from notation to avoid cross-file symbol conflicts

## harvest:const skipped
- Astronomical constants (M_sun, pc, AU): need an IAU source not yet in refs/
- Radiation constant a = 4 sigma/c: not listed in the CODATA tables (only sigma, c1, c1L, c2), so not harvested
- c1L (first radiation constant for spectral radiance): unit cell garbled in the text layer, not needed
- Notation symbols c and G: already declared in eos.notation-structure with different units, not redeclared (see notes)
- CGS conversions: the source states no CGS factors, so all values are SI as printed and CGS is left to the code
- Compton wavelength, molar gas constant, Wien constants, magnetic moments and the other Table XXXI rows: not needed by the PRD

## harvest:obs skipped
- Walter 2010 position angle 100.3+/-0.3 deg (in the statement of the proper-motion entry only, no separate entry)
- Walter 2010 combined WFPC2+HRC parallax 8.33+/-0.6 mas (120 +8/-10 pc)
- Walter 2010 pre-allowance distance 123 +11/-13 pc (recorded in the distance entry's statement only)
- Kaplan 2011 excess at 1500 A (7.1+/-0.4) and extrapolated STMAG for F606W and F140LP (statement text only, no separate entry)
- Kaplan 2011 measured HST/STIS/WFPC2 magnitudes of RX J1856.5-3754: the paper does not tabulate them (Table 1 covers only five other sources; RX J1856 uses van Kerkwijk & Kulkarni 2001 data), so there are no measured per-band magnitudes except Walter 2010 F475W
- Kaplan 2011 ACS zero points (26.67, 25.75, 20.316): Table 1 footnote, not needed for RX J1856 entries
- Burwitz 2003 other-instrument rows of Table 2 (ROSAT, XMM pn/MOS2/RGS), Lbol 4.1e31 erg/s as separate entry, optical-component limits kT_opt<33.6 eV and R_opt>16.3 km (superscript text garbled), one-component R~12.3 km and 13.3 km
- Ho 2007 Table 1 B=3e12 G column (NH 1.18(2), T_inf 4.34(3), zg 0.27(2)) and blackbody column (NH 0.91(1), T_inf 7.36(1), R_inf 5.0): in notes only
- Ho 2007 Table 2 per-latitude B/Teff/Theta_B values: summarised in the assumptions entry's notes only
- van Kerkwijk & Kaplan 2008 nu, nu-dot (-5.98(14)e-16 Hz/s), Edot 3.3e30 erg/s, t0, TOA rms as separate entries
- Spin-down-derived pulsed fraction from Tiengo & Mereghetti 2007 (~1%): only quoted second-hand by Kaplan 2011
- Not stated in any of the five sources: surface temperature map, numerical optical magnitudes in V/B/R bands for RX J1856, any measured B field independent of timing/atmosphere models

## review-r1:colour_srgb
- Mechanical check, run first: `161 entries: {'proposed': 161}; 0 errors, 0 warnings` (exit 0).
- `refs/stokes_1996_srgb.txt:22` explicitly identifies an obsolete proposal with a subsequently corrected numerical error. These verdicts establish historical transcription, not equivalence to the IEC standard.
- GIFs `srgb09–22` were visually inspected, with magnification where needed. Record `image_verified = true` for accepted poor-OCR equations before blessing. Table 0.2 was checked against clean original HTML; its entries can use `ocr = "good"`.
- `colour.toml` duplicates the XYZ notation, matrices, transfer functions, white point and primaries. Its primed inverse conditions conflict with this file's literal unprimed conditions. Consolidate entries and distinguish source transcription from derived corrections.
- Namespace chromaticity `x,y,z` against pixel coordinates and atmospheric formulas. Explicitly map `wavelength ↔ λ` and `xbar,ybar,zbar ↔ \bar{x},\bar{y},\bar{z}`.
- Preserve the unresolved D65 discrepancy: `y=0.3291` at `refs/stokes_1996_srgb.txt:621`, versus `0.3290` at 773. Remove the short-wavelength warning's dependency on the dropped historical-origin entry.

Worker notes: (1) colour.srgb-inverse-power was rejected per the reviewer, but its stated reason is factually wrong. srgb21.gif prints G' and B' in the G and B rows, as I confirmed by viewing it. The real defect is that the printed condition (srgb20.gif) is unprimed and the entry doesn't state it. I recorded this dispute in the review note, and the entry can be re-proposed with the unprimed condition stated. (2) I did not apply the reviewer's fix to the Y meaning in colour.notation-stokes-cie-xyz, which would drop "luminance-like". colour.toml's colour.notation-cie-xyz has the identical Y symbol text, and changing only one makes check.py report a symbol-conflict error. Resolve both together. (3) The short-wavelength-error entry's depends_on was changed from cmf-origin (now rejected) to colour.notation-cie-cmf. (4) Accepted ocr=poor entries got image_verified=true. The reviewer said it inspected the GIFs, and the HTML for the primaries and D65 table.

## review-r1:condensed
- Mechanical check, run first: exit 0; `357 entries: {'proposed': 357}; 0 errors, 215 warnings`. Warnings comprise 212 missing-provenance notices elsewhere and the `alpha`, `lambda`, `phi` conflicts between the two GR notation entries.
- Image verification: all 22 `ocr = "poor"` entries were re-extracted and checked against PDF page images. All 22 lack `image_verified`; no equation remained unreadable.
- Polarization conventions: resolve the 2012 prose/Appendix-B contradiction before mapping linear polarization to Stokes or X/O quantities. The opacity entries use j=1,2 for X/O modes, unlike these linear basis indices. Distinguish reflected-component row sums from incident-component column sums.
- Notation and units: explicitly map `E_ph`/`E`, `B_{14}`/`B_14`, and the several field-normal angle names. Keep that angle separate from opacity's `theta_B` (k-B angle) and surface's magnetic colatitude. Specify Gaussian-CGS constants and the intensity energy-unit conversion; do not silently reuse SI CODATA numbers or the rounded EOS `m_u`.
- Reimplementation: specify the quartic variable as n_{z,j} and apply the decay condition to its normal component. The source distinguishes this from scalar refractive index. [2012:725-753,819-830]
- Need: no clearly out-of-scope entry warrants `drop`; these support Milestone 3. Fits and their truncation repairs remain verification targets under rule 7.
- Session records: the read-only filesystem prevented HANDOFF/worklog updates.

## review-r1:const
- The schema and database changed during review. A later run reported `357 entries: {'proposed': 357}; 0 errors, 215 warnings`. All 30 constant entries lack the now-required `provenance` field; this blocks blessing independently of the scientific verdicts above.
- Symbol spelling hides collisions: fine-structure `α` also appears as geometrical `alpha` and photometric `\alpha` (`gr.toml:17,88`; `obs.toml:618`). Preserve explicit source-qualified meanings.
- Complete the alias map for `ℏ`/`hbar`, `α`/`alpha_f`, and `σ`/`σ_SB`/`σ_sb` (`opacity.toml:39,48`; `surface.toml:119`; `rt.toml:100`).
- Constants here use SI; EOS and opacity notation uses CGS. In particular, `e` is in coulombs here and esu at `opacity.toml:37`. No sourced conversion entry currently bridges those conventions.

(Applier note: the 'fix' entries were each checked against refs/tiesinga_2021_codata2018.txt and applied as given. The reject entry is const.notation-fundamental; the drops are const.avogadro and const.neutron-mass. None of the const entries has ocr = "poor". check.py reports 0 errors; the warnings are the missing provenance/ownership fields.)

## review-r1:obs
- **Mechanical check, first run:** exit 1; `318 entries: {'proposed': 318}; 10 errors, 250 warnings`. Errors were the short `opacity.coulomb-log-classical` anchor and undeclared conflicts for `Z`, `e`, `u`, `rho`, `xi`, `b`, `B_12`, `a_m`, `nu`. A subsequent live-tree run returned exit 0: `357 entries: {'proposed': 357}; 0 errors, 215 warnings`. `obs.toml` remained unchanged during this audit.
- **Blessing remains blocked:** all 36 non-notation entries lack required `provenance`. Scientific accepts above do not waive that schema requirement. Distinguish observations, conventions, inferred quantities and external model outputs; supply `reimplement` or `drop_reason` wherever required.
- **Notation escapes the checker:** `alpha`/`α`/`\alpha` and `nu`/`ν`/`\nu` conceal collisions. Explicitly map the radius/temperature spelling variants, `N_H^{fit}`/`N_H`, and `P_{1856}`/`P_spin`, preserving units and model context.
- **Dependencies and independence:** measured spin period does not depend on dipole braking. Split the rejected composite inference entry. Walter's two means reuse the same observations; Kaplan's quoted X-ray parameters reuse Burwitz. Ho's fitted parameters must retain their shared atmosphere assumptions.
- **R5 remains under-specified:** the Walter magnitude lacks a confirmed system; Kaplan's 27.13 is a prediction and Table 5 contains fitted summaries. Add measured optical/UV band targets and checksummed calibration inputs. Preserve Ho's 0.12-mag extinction separately from Kaplan's 0.05-mag choice.

## review-r2:colour_srgb
- First mechanical check, exit 0: `357 entries: {'blessed': 32, 'rejected': 11, 'proposed': 314}; 0 errors, 532 warnings`.
- `truth/entries/colour.toml:291` duplicates the unsupported "luminance-like" description. Consolidate XYZ notation; the accepted colorimetric-basis entry still depends on the rejected notation entry.
- Record `image_verified = true` for the three reviewed poor-OCR entries. This review inspected the local HTML-linked GIFs under `refs/data/stokes_1996_srgb_img/`; no corresponding local PDF exists.
- The four reviewed non-notation entries lack `provenance` and `ownership`. Scientific acceptance does not resolve those schema omissions.
- `refs/stokes_1996_srgb.txt:22` identifies an obsolete proposal with a subsequently corrected numerical error. These acceptances establish historical transcription, not IEC equivalence; preserve the printed-versus-interpreted decode distinction.

## review-r2:condensed
- **Mechanical check, run first:** exit 0; `357 entries: {'blessed': 220, 'rejected': 46, 'proposed': 91}; 0 errors, 524 warnings`. Warnings: 208 missing provenance, 313 missing ownership, three GR symbol conflicts (`alpha`, `lambda`, `phi`).
- **Image verification:** all 21 reviewed `ocr = "poor"` entries were independently checked against page images; fresh page extractions matched the local text. All were readable.
- **Polarization:** resolve the §2.2.2/Appendix-B basis contradiction before Stokes or X/O mapping. [2012:215–216](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.txt:215), [2012:841–858](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.txt:841).
- **Notation/units:** explicitly connect the remaining photon-energy, field-scaling and field-normal-angle aliases. Specify Gaussian-CGS constants and per-energy intensity conversions; the CODATA entries supply SI values.
- **Reimplementation:** specify the quartic variable as \(n_{z,j}\) and apply decay selection to the normal component; `Im n<=0` remains ambiguous. [2012:725–753,819–830](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.txt:725).
- **Need:** all reviewed claims support Milestone 3; published fits remain verification targets.
- **Session records:** the read-only filesystem prevented HANDOFF/worklog updates.

## review-r2:const
- The three reviewed constants still lack `provenance` and `ownership`; scientific acceptance does not resolve these schema warnings.
- Alias mapping remains incomplete for `ℏ`/`hbar`, `α`/`alpha_f`, and `σ`/`σ_SB`/`σ_sb`. Mass and radiative-constant notation also spans SI and CGS; conversions need explicit provenance.
- Fine-structure `α` visually collides with geometrical `alpha` (`gr.toml:17,88`) and photometric `\alpha` (`obs.toml:618`). Declare source-qualified meanings and the adopted project symbols.

## review-r2:obs
- **Mechanical check, first run:** exit 0; `357 entries: {'blessed': 126, 'rejected': 23, 'proposed': 208}; 0 errors, 532 warnings`. Warning categories: 212 missing provenance, 317 missing ownership, three GR symbol conflicts (`alpha`, `lambda`, `phi`).
- **Metadata remains incomplete:** all six reviewed non-notation entries lack `provenance` and `ownership`. Scientific acceptance does not waive those schema requirements.
- **Remaining synonyms:** explicitly reconcile `N_H^{fit}`/`N_H`, `P_{spin}`/`P_spin`, and apparent-radius/temperature spellings across domains. Preserve km/cm conversions and emission-region versus whole-star distinctions.
- **Dependency error:** the measured spin period still depends on the rejected dipole-braking entry. Remove that dependency; the timing measurement does not require the braking inference.
- All nine entries serve Milestones 1-3.