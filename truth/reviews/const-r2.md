`python3 -I truth/check.py` exited 0. Initial output summary:

```text
357 entries: {'blessed': 57, 'rejected': 14, 'proposed': 286}; 0 errors, 532 warnings
```

All source line numbers below refer to `refs/tiesinga_2021_codata2018.txt`.

| id | verdict | reason |
|---|---|---|
| const.hbar | accept | Anchor L3453 is on Table XXXI p. 45 and confirms `ℏ = 1.054 571 817… × 10^-34 J s`, exact. The unsupported defining relation was removed. |
| const.atomic-mass-constant | accept | Anchor L3735 confirms `1.660 539 066 60(50) × 10^-27 kg` and relative uncertainty `3.0 × 10^-10`. L3737 gives `m_u = 2hcR_∞/(α²c²A_r(e))`; footnote f, L3786–3787, gives `m_u = m(12C)/12 = 1 u`. |
| const.stefan-boltzmann | accept | Anchor occurs at the specified Table XXXI row, L3762: `σ = (π²/60)k⁴/(ℏ³c²)`, `5.670 374 419… × 10^-8 W m^-2 K^-4`, exact. The unsupported radiation-constant claim was removed. |
| const.notation-particles | accept | Anchor L3497 confirms `m_e`; L3586, L3590 and L3616 confirm `m_p`, `m_p/m_e` and `m_n`, with matching units. L3786–3787 confirms `u = m_u`. The expanded locator covers every symbol. |
| const.notation-atomic | accept | Anchor L3481 is correctly located. L3469 supports `μ_B`; L3481–3488 supports `α`, `cR_∞`, `hcR_∞`, `R_∞`, `a_0`, `E_h`; L3525–3526 supports `r_e = α²a_0` and `σ_e = (8π/3)r_e²`. L3441–3442 identifies dimensionless `u_r`. The expanded locator includes p. 46. |

- The three reviewed constants still lack `provenance` and `ownership`; scientific acceptance does not resolve these schema warnings.
- Alias mapping remains incomplete for `ℏ`/`hbar`, `α`/`alpha_f`, and `σ`/`σ_SB`/`σ_sb`. Mass and radiative-constant notation also spans SI and CGS; conversions need explicit provenance.
- Fine-structure `α` visually collides with geometrical `alpha` (`gr.toml:17,88`) and photometric `\alpha` (`obs.toml:618`). Declare source-qualified meanings and the adopted project symbols.