| id | verdict | reason |
|---|---|---|
| colour.notation-stokes-cie-xyz | fix | Replace the normalisation sentence with “The source does not specify XYZ normalisation.” Neither arbitrary input scaling nor “luminance-like” is established at `refs/stokes_1996_srgb.txt:786–792`. |
| colour.notation-srgb-linear | accept | Symbols match GIF `srgb09`; the conditional clipping qualification matches `refs/stokes_1996_srgb.txt:786–802`. |
| colour.notation-srgb-nonlinear | accept | Primed channel symbols match GIFs `srgb11/13` and their introduction at `refs/stokes_1996_srgb.txt:805–819`. |
| colour.notation-srgb-digital | accept | WDC/KDC meanings and eight-bit channel interpretation match `refs/stokes_1996_srgb.txt:830–860` and GIFs `srgb14–17`. |
| colour.notation-chromaticity | fix | Table values match `refs/stokes_1996_srgb.txt:755–779`; all four column sums equal 1.0000. Add explicit scope: “Here x,y,z denote chromaticities, distinct from image and atmospheric coordinates.” |
| colour.notation-cie-cmf | accept | Names, column order and wavelength units match `refs/data/cie1931_2deg_cvrl_README.txt:2,6–12`. |
| colour.srgb-colorimetric-basis | fix | Correct locator to “Part 2, ‘A Single RGB Standard Color Space’; ‘Colorimetric definitions and digital encodings’.” The anchor precedes the claimed subsection: `refs/stokes_1996_srgb.txt:588–601`; viewing dependence is at 742–747. |
| colour.srgb-bt709-relation | drop | Historical standards relationship is unnecessary for Milestones 1–3; retain the actual primaries separately. `refs/stokes_1996_srgb.txt:964–969` supports the prose relationship but prints no supplied gamma-containing formula. |
| colour.srgb-stokes-xyz-to-linear-matrix | accept | All nine coefficients, signs and channel ordering match visually inspected GIF `srgb09`; anchor and equation label correspond at `refs/stokes_1996_srgb.txt:786–792`. |
| colour.srgb-stokes-linear-to-xyz-matrix | accept | All nine coefficients match visually inspected GIF `srgb22`, including `0.0193`. The anchor introduces the reverse-equation group at `refs/stokes_1996_srgb.txt:857–881`. |
| colour.srgb-encode-threshold | fix | GIFs `srgb10/12` confirm `0.00304` and both inequalities, corresponding to `refs/stokes_1996_srgb.txt:805–819`. Replace the unsupported IEC comparison with “This is the 1996 proposal threshold.” |
| colour.srgb-stokes-encode-linear | accept | GIFs `srgb10/11` confirm `≤0.00304`, factor `12.92` and corresponding channel inputs; `refs/stokes_1996_srgb.txt:805–812`. |
| colour.srgb-stokes-encode-power | accept | GIFs `srgb12/13` confirm `>0.00304`, factor `1.055`, exponent `1.0/2.4`, subtraction `0.055` and corresponding channel inputs; `refs/stokes_1996_srgb.txt:805–819`. |
| colour.srgb-gamma22-fit | drop | Historical motivation adds no required transform or validation target for Milestones 1–3. `refs/stokes_1996_srgb.txt:823–827` contains the prose fit claim, without the supplied power-law formula. |
| colour.srgb-digital-scaling | accept | GIF `srgb14` confirms `(WDC−KDC)R′+KDC` independently for each channel; definitions and label match `refs/stokes_1996_srgb.txt:830–836`. |
| colour.srgb-digital-8bit | accept | Counts, eight-bit qualification and all channel equations match GIF `srgb15` and `refs/stokes_1996_srgb.txt:840–846`. |
| colour.srgb-digital-8bit-simplified | accept | GIF `srgb16` prints `255.0` times each corresponding encoded channel; anchor matches `refs/stokes_1996_srgb.txt:850–853`. No rounding rule is asserted. |
| colour.srgb-digital-decode | accept | GIF `srgb17` confirms division by `255.0` for each corresponding channel; `refs/stokes_1996_srgb.txt:857–860`. |
| colour.srgb-decode-threshold | fix | GIFs `srgb18/20` confirm `0.03928` and unprimed conditions, corresponding to `refs/stokes_1996_srgb.txt:857–874`. Correct the note to `12.92×0.00304=0.0392768≈0.03928`; explicitly retain the unresolved source-symbol error. |
| colour.srgb-inverse-linear | fix | State the printed condition explicitly: `R_{sRGB},G_{sRGB},B_{sRGB}\le0.03928` **unprimed**, followed by the stated per-channel division. GIFs `srgb18/19`, `refs/stokes_1996_srgb.txt:864–867`; interpreting the condition as encoded input requires a labelled correction. |
| colour.srgb-inverse-power | reject | GIF `srgb21` prints `R′` in the G and B rows too; “same for G,B” silently repairs the source. GIF `srgb20` also uses unprimed conditions. `refs/stokes_1996_srgb.txt:871–874` does not support the claimed complete per-channel transcription. |
| colour.srgb-clipping | accept | The encoder limitation and conditional clipping qualification faithfully match `refs/stokes_1996_srgb.txt:796–802`. |
| colour.srgb-stokes-d65-white-table02 | accept | Triple matches `refs/stokes_1996_srgb.txt:764–779` and original HTML table cells. The differing Table 0.1 value is explicitly disclosed. |
| colour.srgb-d65-white-table01 | accept | Exact pair appears at `refs/stokes_1996_srgb.txt:620–621`; attribution to D65 is repeated at 664–665. |
| colour.srgb-red-primary | accept | Triple matches `refs/stokes_1996_srgb.txt:758–779` and the original HTML Red column; introductory anchor is relevant at 750–752. |
| colour.srgb-green-primary | accept | Original HTML separates the concatenated text cells and confirms `(0.3000,0.6000,0.1000)`; `refs/stokes_1996_srgb.txt:750–779`. |
| colour.srgb-blue-primary | accept | Original HTML separates the concatenated text cells and confirms `(0.1500,0.0600,0.7900)`; `refs/stokes_1996_srgb.txt:750–779`. |
| colour.srgb-viewing-environment | accept | Every stated table value and unit matches `refs/stokes_1996_srgb.txt:608–649`; the recommendation in the notes matches 690–703. |
| colour.cie1931-2deg-cmf-file | accept | `refs/data/cie1931_2deg_cvrl_README.txt:2` supplies spacing and columns. CSV inspection confirms 471 rows, endpoints 360/830 nm, unit steps and twelve decimal places; stated SHA-256 matches `refs/CHECKSUMS:56`. |
| colour.cie1931-2deg-cmf-origin | drop | Historical reconstruction is unnecessary for Milestones 1–3. `refs/data/cie1931_2deg_cvrl_README.txt:16–24` supports the prose, but does not print the supplied barred-function equation or coefficient symbols. |
| colour.cie1931-2deg-cmf-short-wavelength-error | accept | Explicitly attributed warning faithfully matches `refs/data/cie1931_2deg_cvrl_README.txt:24–30`, including the 460 nm qualification and proposed modifications. |

- Mechanical check, run first: `161 entries: {'proposed': 161}; 0 errors, 0 warnings` (exit 0).
- `refs/stokes_1996_srgb.txt:22` explicitly identifies an obsolete proposal with a subsequently corrected numerical error. These verdicts establish historical transcription, not equivalence to the IEC standard.
- GIFs `srgb09–22` were visually inspected, with magnification where needed. Record `image_verified = true` for accepted poor-OCR equations before blessing. Table 0.2 was checked against clean original HTML; its entries can use `ocr = "good"`.
- `colour.toml` duplicates the XYZ notation, matrices, transfer functions, white point and primaries. Its primed inverse conditions conflict with this file’s literal unprimed conditions. Consolidate entries and distinguish source transcription from derived corrections.
- Namespace chromaticity `x,y,z` against pixel coordinates and atmospheric formulas. Explicitly map `wavelength ↔ λ` and `xbar,ybar,zbar ↔ \bar{x},\bar{y},\bar{z}`.
- Preserve the unresolved D65 discrepancy: `y=0.3291` at `refs/stokes_1996_srgb.txt:621`, versus `0.3290` at 773. Remove the short-wavelength warning’s dependency on the dropped historical-origin entry.