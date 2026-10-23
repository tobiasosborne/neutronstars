Initial `python3 -I truth/check.py` output:

```text
161 entries: {'proposed': 161}; 0 errors, 0 warnings
```

All `L…` references below identify lines in [refs/tiesinga_2021_codata2018.txt](/home/tobias/Projects/neutronstars/refs/tiesinga_2021_codata2018.txt). Verdicts concern scientific support; subsequent schema blockers are listed afterward.

| id | verdict | reason |
|---|---|---|
| const.speed-of-light | accept | L3444 confirms `299 792 458 m s^-1`, exact, at the stated table/page; defining status is supported by L119–123. |
| const.newtonian-gravitation | accept | L3449 confirms `6.674 30(15) × 10^-11 m^3 kg^-1 s^-2` and relative uncertainty `2.2 × 10^-5`. |
| const.planck | accept | L3451 confirms `6.626 070 15 × 10^-34 J Hz^-1`, exact; the SI revision is supported by L96, L110–126. |
| const.hbar | fix | L3453 confirms the exact, truncated value `1.054 571 817… × 10^-34 J s`. Remove “(h divided by 2 pi)”: that defining relation is not stated at the locator. |
| const.elementary-charge | accept | L3462 confirms `1.602 176 634 × 10^-19 C`, exact. This supports the SI value. |
| const.boltzmann | accept | L3731 confirms symbol `k`, value `1.380 649 × 10^-23 J K^-1`, exact. |
| const.avogadro | drop | L3730 confirms `6.022 140 76 × 10^23 mol^-1`, exact; no molar calculation is required by PRD Milestones 1–3. |
| const.electron-mass | accept | L3497 confirms `9.109 383 7015(28) × 10^-31 kg` and relative uncertainty `3.0 × 10^-10`. |
| const.electron-mass-energy-mev | accept | L3499–3500 explicitly associate `0.510 998 950 00(15) MeV` with `m_e c^2`; the accompanying joule value also matches. |
| const.proton-mass | accept | L3586 confirms `1.672 621 923 69(51) × 10^-27 kg` and relative uncertainty `3.1 × 10^-10`. |
| const.neutron-mass | drop | L3616 confirms `1.674 927 498 04(95) × 10^-27 kg`; free-neutron microphysics is outside the exterior-spacetime and surface-model work in Milestones 1–3. |
| const.proton-electron-mass-ratio | accept | L3590 confirms `m_p/m_e = 1836.152 673 43(11)`, dimensionless, with relative uncertainty `6.0 × 10^-11`. |
| const.atomic-mass-constant | fix | Replace “footnote f gives” with “the table row gives”: the formula is at L3737. Footnote f, L3786–3788, establishes `m_u = m(12C)/12 = 1 u`; L3735 confirms the numeric value and uncertainty. |
| const.unified-atomic-mass-unit | accept | Anchor is correctly at L3893, Table XXXIV, p. 51: `1.660 539 066 60(50) × 10^-27 kg`. Equality with `m_u` is explicit at L3786–3787. |
| const.electron-volt | accept | L3891 confirms `(e/C) J`, `1.602 176 634 × 10^-19 J`, exact. |
| const.stefan-boltzmann | fix | L3762 confirms the complete formula and `5.670 374 419… × 10^-8 W m^-2 K^-4`, exact. Delete “The source lists no radiation constant a = 4 sigma/c.” The embedded relation for `a` is unsupported at this locator. |
| const.first-radiation-constant | accept | L3765 confirms `2πhc²`, `3.741 771 852… × 10^-16`, exact. PDF p. 49 visually confirms `[W m²]`; footnote h is at L3792. |
| const.second-radiation-constant | accept | L3766 confirms `c_2 = hc/k = 1.438 776 877… × 10^-2 m K`, exact; the unit qualification is at L3774–3776. |
| const.vacuum-permeability | accept | L3445–3446 confirm `4παℏ/(e²c)`, `1.256 637 062 12(19) × 10^-6 N A^-2`, and departure from the former exact value. |
| const.vacuum-permittivity | accept | L3447 confirms `1/(μ_0 c²) = 8.854 187 8128(13) × 10^-12 F m^-1`. |
| const.fine-structure | accept | L3481 confirms the SI expression `e²/(4πε_0ℏc)`, value `7.297 352 5693(11) × 10^-3`, and relative uncertainty `1.5 × 10^-10`. |
| const.classical-electron-radius | accept | L3525 confirms `α²a_0 = 2.817 940 3262(13) × 10^-15 m`. |
| const.thomson-cross-section | accept | L3526 confirms source symbol `σ_e`, factor `(8π/3)r_e²`, and `6.652 458 7321(60) × 10^-29 m²`. |
| const.bohr-radius | accept | L3487 confirms both algebraic expressions and `5.291 772 109 03(80) × 10^-11 m`; denominators were checked against PDF p. 45. |
| const.rydberg-constant | accept | L3486 confirms `10 973 731.568 160(21) m^-1`; L3774 specifies cycles or periods per metre. |
| const.rydberg-frequency | accept | L3483 confirms both expressions, `α²m_e c²/(2h) = E_h/(2h)`, and `3.289 841 960 2508(64) × 10^15 Hz`. |
| const.rydberg-energy-joule | accept | L3484 confirms `hcR_∞ = 2.179 872 361 1035(42) × 10^-18 J`. |
| const.rydberg-energy-ev | accept | L3484–3485 associate `13.605 693 122 994(26) eV` with the same `hcR_∞` energy equivalent. |
| const.hartree-energy | accept | L3488 confirms all three defining expressions and `4.359 744 722 2071(85) × 10^-18 J`; PDF p. 45 confirms denominator grouping. |
| const.bohr-magneton | accept | L3469 confirms `eℏ/(2m_e) = 9.274 010 0783(28) × 10^-24 J T^-1`. |
| const.notation-fundamental | reject | Anchor L3449 identifies `G`, which the symbol list excludes. It does not substantiate that list; several listed quantities instead occur on pp. 49 and 51. `h_P` is project notation, whereas the statement claims source notation. Split and re-anchor. |
| const.notation-particles | fix | Replace locator with `"Table XXXI, pp. 025010-45, 025010-47, 025010-49 (footnote f)"`. Electron mass is at L3497; proton, neutron and ratio at L3586, L3616, L3590; `u = m_u` at L3786–3787. |
| const.notation-atomic | fix | Replace locator with `"Table XXXI, pp. 025010-45–025010-46"`: `r_e` and `σ_e` occur at L3525–3526 on p. 46, outside the current locator. Other listed symbols match p. 45. |

- The schema and database changed during review. A later run reported `357 entries: {'proposed': 357}; 0 errors, 215 warnings`. All 30 constant entries lack the now-required `provenance` field; this blocks blessing independently of the scientific verdicts above.
- Symbol spelling hides collisions: fine-structure `α` also appears as geometrical `alpha` and photometric `\alpha` (`gr.toml:17,88`; `obs.toml:618`). Preserve explicit source-qualified meanings.
- Complete the alias map for `ℏ`/`hbar`, `α`/`alpha_f`, and `σ`/`σ_SB`/`σ_sb` (`opacity.toml:39,48`; `surface.toml:119`; `rt.toml:100`).
- Constants here use SI; EOS and opacity notation uses CGS. In particular, `e` is in coulombs here and esu at `opacity.toml:37`. No sourced conversion entry currently bridges those conventions.