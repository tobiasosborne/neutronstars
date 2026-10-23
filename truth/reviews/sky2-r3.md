| id | verdict | reason |
|---|---|---|
| obs.mellinger-photometric-calibration | reject | The paraphrase faithfully combines `refs/mellinger_2009_panorama.txt:285–288`; both calibration equations match lines 290–302, and TeX anchors match byte-for-byte. However, the round-2 locator correction remains unapplied: add **Sects. 2.1.3 and 3.5** for filters and TIFF processing (lines 93–95, 342–346). **Surviving statement:** typically 100–120 YBS matches per field/channel, with calibration given by Eqs. (4)–(5). |
| obs.rxj1856-epoch-propagation-estimate | accept | Locator correction to **(D.11)** is complete: `refs/butkevich_lindegren_2014_propagation.txt:1364–1377`. TeX anchor matches at line 2291; the reduced estimate follows (D.9), text lines 1345–1354, using Walter Table 3 means (`refs/walter_2010_rxj1856_distance.txt:862–874`). Approximation and illustrative timespan remain explicit. |

- Mechanical check, run first, exit 0: `493 entries: {'blessed': 392, 'rejected': 96, 'proposed': 5}; 0 errors, 0 warnings`.
- Independent recomputation: expected `Δθ = 4.36884266e-5 mas`; computed `4.368842659710374e-5 mas`; absolute error `2.896e-15 mas`; tolerance `5e-12 mas`. Equation: BL (D.9), cited above.
- No new cross-entry notation issues from these two corrections.