# Notation merge report, 2026-10-23 (lane notation-merge)

Input: `truth/harvest/notation_additions_<lane>.md` for geodesic, camera, const, radiance, field_sky, sky, sky2 and display2. They were read twice. The second read, at 13:10, picked up the revisions that the camera, sky and display2 lanes made after review (camera bridges, the gaia-photometry split, the psf-scoped display symbols), and that read is the one merged.
Output: 102 new rows in `truth/NOTATION_DECISIONS.md` (Dictionary), and a section "## Additions 2026-10-23 (Milestone-1 harvest)" that records the decisions, the collision table, the alias bridges and the new doubts.
Edits: only `truth/NOTATION_DECISIONS.md` and this file. No `truth/entries/*.toml` was touched.

## Counts

| lane | proposed lines | new rows, as proposed | new rows, renamed | merged into another lane's row | bridge to an existing row / reuse | rejected |
|---|---|---|---|---|---|---|
| geodesic | 16 rows + 5 alias changes | 15 | 1 (`\Theta_K(\theta)` -> `\Theta_K(\theta_{Sch})`) | 0 | 5 alias/note changes | 0 |
| camera | 24 rows + 6 bridge remarks | 24 | 0 | 0 | 6 remarks entered in the collision table | 0 |
| const | 6 lines (q_esu and q_C on one line) | 3 (`\mathrm{au}`, `\mathrm{pc}`, `(GM)^{N}_\odot`) | 0 | 0 | 0 | 3 lines (`R^{N}_\odot`, `L^{N}_\odot`, `q_{esu}`/`q_{C}`): entries rejected |
| radiance | 8 | 2 (`B_\lambda`, `K_{cd}`) | 0 | 0 | 5 (two sigma lines, `k_B`, `D`, `R_{em}^\infty`) | 1 (call form `R^\infty(123 pc)`) |
| field_sky | 10 | 4 (`a_{dip}`, `\theta_{surf}`, `\phi_{surf}`, `\mu_{mag}`) | 0 | 2 (`\alpha_{RA}`, `\delta_{Dec}`, merged with sky; meaning made generic) | 3 (`\theta_B` <- chi, `\theta_{Bn}` <- theta, `x_g` <- x) | 1 (`\alpha_B`, `i`) |
| sky | 17 | 14 (`\alpha_{RA}`, `\delta_{Dec}` with the generic meaning) | 2 (`\langle f_\lambda\rangle`, `\langle f_\nu\rangle` -> `\langle F_\lambda\rangle`, `\langle F_\nu\rangle`) | 1 (`A_{u}` -> the const row `\mathrm{au}`) | 0 | 0 |
| sky2 | 22 | 22 | 0 | 0 | 0 | 0 |
| display2 (revised) | 17 | 12 | 3 rows from 2 lines (`x_1` -> `x_{1,psf}`; `K_\lambda`, `A_\lambda` -> `K_{psf}`, `A_{psf}`) | 0 | 3 (`L_{sc}`, `L_d`, `\lambda` reused) | 0 |
| **total** | **120 lines** (+5 alias changes, +6 bridge remarks) | **96** | **6** | **3** | **16** | **5 lines** |

New Dictionary rows: 96 + 6 = **102**. The dictionary went from 369 to 471 distinct `sym`, with no duplicates and no literal clash with an existing row or with the retired bare symbols.
Symbols renamed relative to the proposals: 7. These are `\Theta_K(\theta)`, `A_{u}`, `\langle f_\lambda\rangle`, `\langle f_\nu\rangle`, `x_1`, `K_\lambda` and `A_\lambda`. The first draft of display2 also had `I_0`, `\theta_1`, `L_{in}` and `L` -> `L_{abs}`; the lane itself replaced those with `I_{peak}`, `\theta_{1,psf}`, `L_{sc}` and `Y_{pre}` before the second read.

Rejected, with reasons:
- `R^{N}_\odot`, `L^{N}_\odot` (const): const.r-sun-nominal and const.l-sun-nominal are status rejected. The gpt "drop" verdict says no Milestone 1-3 step uses them.
- `q_{esu}`, `q_{C}` (const): const.esu-from-coulomb is status rejected (gpt r1 and r2). NIST prints only the rounded 3.335 641e-10 C per statC, so the exact 10 c conversion is not identified with Gaussian charge.
- `R^\infty(123\,\mathrm{pc})` (radiance): Ho07's radius is the emission-region radius `R_{em}^\infty`, as the existing collision row `R_inf` and the lane's own later row both say. The call form is `R_{em}^\infty(123\,\mathrm{pc})`.
- `\alpha_B`, `i` (field_sky): these are informal names for the existing `\Psi_{pole}` and `\theta_{los}`, and bare `i` is a source spelling of the mode index `j`.

## Collision table (summary; the full table is in NOTATION_DECISIONS.md, "Additions 2026-10-23")

| kind | sym / source spelling | lanes | decision | one-line reason |
|---|---|---|---|---|
| (a) same sym, different meaning | `\alpha_{RA}`, `\delta_{Dec}` | sky (Gaia ICRS, J2016.0) vs field_sky (Walter+2010, epoch not stated) | one row with a generic meaning | Frame and epoch are attributes of a datum, not of the symbol. |
| (b) two syms, one quantity | `A_{u}` (sky) vs `\mathrm{au}` (const) | sky, const | `\mathrm{au}` | Unit-style spelling like `\mathrm{eV}`, `\mathrm{pc}`, `\mathrm{Ry}`; IAU's own symbol. |
| (b) | `\tilde E`, `\tilde L` vs the lead's `E_{geo}`, `L_{geo}` | geodesic | `\tilde E`, `\tilde L` | `geo` means geometric units (`u_{geo}`); the family matches the reserved `\tilde E_{cam}`, `\tilde L_{cam}`. |
| (b) | `\tau_{prop}` vs the reserved `\tau_{cam}` | geodesic, pipeline table | both | Generic proper time, and its camera value. |
| (b) | `r_{peri}` vs PRD `r_p`; `r_{cam}` vs pipeline `r_c` | geodesic | `r_{peri}`, `r_{cam}` | Subscript p means proton/polar; `r_{cam}` is reserved. |
| (b) | `f_{cam}` (pixels) vs `f_{foc}` (cm) | camera, display2 | both kept | Same lens in different units: f_cam = f_foc / pixel pitch (likewise r_pix = r_fp / pixel pitch). |
| (b) | `u^\mu_{obs}` (Pech83, blessed) vs label `st` (camera) | camera | keep `st` | `\nu_{obs}` would read as the camera frequency; unifying them would need a blessed-entry edit (doubt). |
| (b) | `\langle f_\lambda\rangle` vs the existing `F_\lambda` | sky | `\langle F_\lambda\rangle`, `\langle F_\nu\rangle` | Flux density is F in this dictionary; `f_\nu` is the Eddington factor. The sky-r1 reviewer suggested the same. |
| (b) | `L_{sc}` (existing Reinhard row) vs display `L_{in}`/`L` | display2 | reuse `L_{sc}`, new `Y_{pre}` | The lane's revision; the Reinhard key form is a special case of the exposure convention. |
| (c) existing row | `I_0` (Bel02) vs Airy central intensity | display2 | `I_{peak}` | `I_0` is taken. |
| (c) | `J_{l}` source spelling `J_1` vs the Bessel function | display2 | `\mathrm{J}_1` | Upright special function. |
| (c) | `\theta` (RT ray angle) vs DA09, James, MG09, PY01, Airy, BL14 thetas | geodesic, camera, field_sky, display2, sky2 | `\theta_{Sch}`, `\theta_{cam}`, `e^{(a)}`, `\theta_{Bn}`, `\theta_{psf}`, `\Delta\theta_{lt}`; `\Theta_K(\theta_{Sch})` | Bare theta is the RT ray angle. |
| (c) | `x` (CIE) vs PY01, WS83, Airy x; `x_1` | field_sky, display2 | `x_g`, `x_{psf}`, `x_{1,psf}`; WS83 x in latex only | Bare x is chromaticity; x_1 takes the psf scope of `\theta_{1,psf}`. |
| (c) | `A_V` (extinction) vs BL14 `A_V`, Airy `A_\lambda` | sky2, display2 | `A_{vel}`, `A_{psf}` (and `K_{psf}`) | A_lambda reads as extinction, and subscript lambda means per unit wavelength (`F_\lambda`, `B_\lambda`, `I_\lambda`). |
| (c) | `\varpi_\nu` albedo vs parallax | sky | `\varpi_{plx}` | Keeps the glyph off the albedo rows. |
| (c) | `G`, `P`, `S`, `n_p`, `\mu`, `b`, `r`, `p`, `q`, `u`, `tau`, `D`, `g`, `a`, `\beta`, `\gamma`, `\Lambda`, `E`, `R`, `L_z` | all lanes | subscripted or decorated as in the rows | Each bare letter is taken or retired; see the dictionary collision table. |
| (c) retired | `L`, `u`, `n`, `S`, `\xi` | display2, sky2, camera, sky | `Y_{pre}`/`L_{sc}`, `\vec{u}_{bary}`/`u^\mu_{cam}`, `\hat v_{cam}`, `S_{Gaia}`; MG09 `\xi` latex only | Retired list. |
| source-level | `nu_0` | camera (DA09 observed) vs Bel02 (emitter) | DA09 -> `\nu_{cam}`, Bel02 -> `\nu` | The same source spelling has opposite frames. |

## Entry edits required (for the apply passes)

Notation entries (check.py's "one sym, one row" rule). The two `\langle F\rangle` renames were applied by the apply-sky-r2 lane during this merge, following the sky r2 review's identical decision. The only sym still matching no row is `A_{u}`:

| file | entry id | field | old | new |
|---|---|---|---|---|
| truth/entries/sky.toml | obs.notation-gaia-sky | symbols[].sym | `A_{u}` | `\mathrm{au}` |
| truth/entries/sky.toml | obs.notation-gaia-sky | symbols[`\mathrm{au}`].meaning | `astronomical unit (as used in mu_r = v_r varpi / A_u)` | `astronomical unit, exactly 149 597 870 700 m (IAU 2012 B2)` (keep `sym_source = "A_u"`, units `cm`; consider `depends_on` const.au) |
| truth/entries/sky.toml | obs.notation-gaia-photometry | symbols[].sym | `\langle f_\lambda\rangle` | `\langle F_\lambda\rangle`: DONE (apply-sky-r2) |
| truth/entries/sky.toml | obs.notation-gaia-photometry | symbols[].sym | `\langle f_\nu\rangle` | `\langle F_\nu\rangle`: DONE (apply-sky-r2) |

Meaning edits, with `sym` unchanged. The dictionary row is generic so that the Walter+2010 datum can share it:

| file | entry id | sym | old meaning | new meaning |
|---|---|---|---|---|
| truth/entries/sky.toml | obs.notation-gaia-sky | `\alpha_{RA}` | `right ascension (ICRS, Gaia-CRF3) of a source at the reference epoch J2016.0` | `right ascension of a source; the reference frame and epoch are attributes of each datum (Gaia DR3: ICRS, Gaia-CRF3, epoch J2016.0; Walter+2010: not stated)` |
| truth/entries/sky.toml | obs.notation-gaia-sky | `\delta_{Dec}` | `declination (ICRS, Gaia-CRF3) of a source at the reference epoch J2016.0` | `declination of a source; the reference frame and epoch are attributes of each datum (Gaia DR3: ICRS, Gaia-CRF3, epoch J2016.0; Walter+2010: not stated)` |

Every sym and meaning in obs.notation-gaia-propagation already matches its row (22 of 22), and so does every one in obs.notation-gaia-photometry (10 of 10). The `G_{Gaia}` row was aligned to the r2 meaning "Gaia G-band magnitude". The `|b|` in `\vec{u}_{bary}` differs only by the markdown pipe escape. All 8 `symbols[]` of rt.airy-intensity-pattern match their rows in sym, meaning, units and adopted entry.

Latex and statement edits in non-notation entries, where the latex is the project's own (a derivation or convention, not a paper spelling):

| file | entry id | field | old | new |
|---|---|---|---|---|
| truth/entries/display_psf.toml | rt.airy-first-zero (proposed) | latex, statement | `x_1` | `x_{1,psf}` |
| truth/entries/display_psf.toml | rt.psf-aperture-convention (blessed, pipeline) | latex, statement | `K_\lambda`, `A_\lambda` | `K_{psf}`, `A_{psf}` (wavelength as argument) |

Recommended, not required:
- truth/entries/surface_field.toml, surface.field-colatitude-from-geometry (blessed pipeline): the latex writes `\chi`, `i`, `\alpha_B` for the sub-observer case. Use `\theta_B`, `\theta_{los}`, `\Psi_{pole}`.
- truth/entries/sky2.toml, obs.notation-gaia-propagation notes: "Reuses ... A_u" should read `\mathrm{au}`. "tau the camera proper time" is wrong: the camera proper time is `\tau_{cam}`, and bare tau is optical depth.

## Retired bare symbols in the new non-notation entries (task 5)

Scope: every non-notation entry in gr_geodesic, camera, const_iau, radiance_display, surface_field, sky, sky2 and display_psf, as of 13:10.

`symbols[].sym` violations: **0**. The only non-notation entry with a `symbols` field is rt.airy-intensity-pattern, and its 8 syms are all decorated.

`latex` hits. These are allowed when they keep the paper's spelling, and all of them do:

| file | entry id | bare symbol | source spelling of | project symbol |
|---|---|---|---|---|
| gr_geodesic.toml | gr.bending-integral-finite-radius | `u` in `(1-2u)u^2 ... du` | Pech83 u = M/r | `u_{geo}`. The same latex mixes in project spellings (`r_{cam}`, `b_{\max}`); if it is rewritten, use `u_{geo}`, `u_{b,geo}`. |
| camera.toml | gr.tetrad-vector-components | `\xi` | MG09 azimuth of a tetrad vector | none (latex only) |
| sky.toml | obs.gaia-band-flux-integral | `S(\lambda)` | Riello passband | `S_{Gaia}` |
| sky.toml | obs.gaia-synthetic-mean-fluxes | `S(\lambda)` | Riello passband | `S_{Gaia}` |
| sky.toml | obs.gaia-vegamag-definition | `S(\lambda)` | Riello passband | `S_{Gaia}` |

Decorated forms of retired letters, reported for completeness (allowed): `\mathbf{u}`, `\mathbf{n}` in gr.camera-four-velocity-static-frame (MG09; project `u^\mu_{cam}`, `\hat v_{cam}`); `n_a`, `n^a` in gr.tetrad-lorentz-transformation (MG09; `\hat v_{cam}`); `\vec{u}` in obs.bl-proper-motion-components, obs.bl-propagation-direction-nolt, obs.bl-propagation-pm-nolt, obs.bl-propagation-direction-rigorous, obs.bl-propagation-pm-rigorous and obs.agis-standard-model-direction (`\vec{\bar u}`) (BL14/AGIS; `\vec{u}_{bary}`); `\mathcal{L}^{\mathrm N}_\odot` in const.l-sun-nominal (rejected entry); upright `\mathrm{H}` (a matrix) in obs.gaia-xp-ecs-spd-sed-units.
Resolved since the first read: colour.display-exposure-convention used bare `L` in its latex. The apply pass rewrote it as `L_{sc} = t_{exp}\,k_{exp}\,Y_{pre}`.

## CODATA sentence (task 6): done as instructed, but the premise is wrong

The sentence was added verbatim after the e[esu] sentence, and no numbers were changed. Re-computing it with 40-digit decimals shows that the 3.3e-11 agreement exists only because the printed hbar is truncated:

| quantity | value | source |
|---|---|---|
| 10 c[m/s] e[C] | 4.803204712570e-10 | e = 1.602176634e-19 C, c = 299792458 m/s (exact) |
| sqrt(alpha hbar c), hbar = 1.054571817e-34 (printed, truncated "...") | 4.8032047124135e-10, gap +3.264e-11 | const.toml hbar entry, alpha = 7.2973525693e-3 |
| sqrt(alpha hbar c), hbar = h/(2 pi) = 1.0545718176461564e-34 | 4.8032047138850e-10, gap **-2.737e-10** | h = 6.62607015e-34 J s (exact) |
| mu_0/(4 pi 1e-7) - 1, halved | 5.444e-10 / 2 = 2.722e-10 | mu_0 = 1.25663706212e-6 N A^-2 (CODATA 2018) |

With the exact hbar (h is exact, so h/2pi is too), the gap is 2.7e-10. That matches the r1 reviewer's figure, and it is explained by mu_0 != 4 pi 1e-7 in the 2019 SI. It is larger than u_r(alpha) = 1.5e-10. So "agree to 3.3e-11 ... inside u_r(alpha)" in the Conventions paragraph does not survive, and the reviewer who called the number incorrect was right; const.esu-from-coulomb's r1 note says the same ("the 4.803204712413e-10 figure truncates hbar"). The added sentence is literally true (the 3.3e-11 gap is an artefact of the rounding), but the paragraph's conclusion is not. **Lead decision needed.** The fix would be to state that 10 c e[C] and sqrt(alpha hbar c) differ by 2.7e-10, from the exact h/2pi and the SI-2019 mu_0, and to drop the "holds to that level" claim. The disputed note on the `e` row and the "Rejected (kept, disputed): `e`" line in "Review r1 outcome" carry the same error.

## Checker

`python3 -I truth/check.py`, last line after the merge:

```
493 entries: {'proposed': 52, 'blessed': 350, 'rejected': 91}; 0 errors, 0 warnings
```
