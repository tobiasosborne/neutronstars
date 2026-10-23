| id | verdict | reason |
|---|---|---|
| surface.py01-gr-dipole-magnitude | accept | R1 corrections verified. Anchor at [PY01.txt:223–224](/home/tobias/Projects/neutronstars/refs/potekhin_yakovlev_2001_cooling.txt:223); TeX:484 matches every factor in eq. (10). Flat limit is algebraically correct; added vacuum assumptions agree with WS83:179–201. |
| surface.py01-gr-dipole-field-angle | accept | Anchor at [PY01.txt:216–219](/home/tobias/Projects/neutronstars/refs/potekhin_yakovlev_2001_cooling.txt:216), equation at TeX:487. Direct field-normal definition at txt:199–200 replaces the former doubt correctly. Component signs support the stated quadrant continuation. |
| surface.py01-gr-dipole-a-of-x | accept | Anchor at [PY01.txt:223](/home/tobias/Projects/neutronstars/refs/potekhin_yakovlev_2001_cooling.txt:223); eq. (11), txt:220–224 and TeX:489–490, matches numerator, denominator, signs and square root. Corrected numerical table and coefficients `1/2, 1/8, 11/160` independently verified. |
| surface.ws83-gr-dipole-components | accept | Fresh extraction remains garbled; independently read PDF pp. 4–5 visually. Both components and asymptotes match [WS83.txt:256–289](/home/tobias/Projects/neutronstars/refs/wasserman_shapiro_1983_gr_dipole.txt:256). Anchor L289 is correctly located. Added conditions and geometrized units match L179–201; dependency removal is correct. |
| surface.field-colatitude-from-geometry | fix | R1 frame, material-point motion and dependency corrections are correct; PRD:33 anchor is valid. Replace “point at the centre of the observed disc” with **“geometric sub-observer surface point”**. Their identification is not generally preserved for the moving camera: this follows from the aberration transformation in [James.txt:1175–1181](/home/tobias/Projects/neutronstars/refs/james_2015_interstellar.txt:1175). The dot-product formula survives unchanged. |

- Mechanical check, run first; exit **0**, exact output: `458 entries: {'proposed': 64, 'blessed': 304, 'rejected': 90}; 0 errors, 0 warnings`.
- Recomputed PY01 eq. (11) using 150-digit arithmetic; 60 significant digits shown:

  ```text
  a(0.1)  = 0.513240103203431705735036817975877537520244399590312644475368
  a(0.01) = 0.501256923822306734143338917832779694156633327281512820404545
  a(1e-4) = 0.500012500687548441316167633775211074283290264209485568684307
  ```

  Against R1’s rounded expectations `0.513240103203431706`, `0.501256923822306734`, `0.500012500687548441`, absolute errors are respectively **2.650×10⁻¹⁹, 1.433×10⁻¹⁹, 3.162×10⁻¹⁹**, below **10⁻¹⁵**. Independently evaluated WS83 component ratios agree within **5.375×10⁻¹⁴²**, tolerance **10⁻¹³⁰**.
- Notation mappings are now explicit, but `a_{dip}`, `\mu_{mag}`, `\theta_{surf}` and `\phi_{surf}` remain absent from canonical notation declarations. Merge the additions; remove their stale claim that PY01’s field-normal identification is merely inferred.
- All five entries serve Milestones 1–3; none warrants dropping.