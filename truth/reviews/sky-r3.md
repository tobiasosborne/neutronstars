| id | verdict | reason |
|---|---|---|
| obs.gaia-synthetic-mean-fluxes | accept | Malformed historical quotation removed. Both equations match `refs/riello_2021_gaia_photometry.txt:1108–1121` and `.tex:628–633`; anchors remain byte-identical. |
| obs.gaia-zero-points-vegamag | accept | Malformed historical quotation removed. Expected/stored values 25.6874, 25.3385, 24.7479, each ±0.0028: errors 0, tolerance 0 (`refs/riello_2021_gaia_photometry.txt:1152`). Printed-ratio qualification retained (lines 1091–1098). |
| obs.gaia-xp-absolute-calibration | accept | Statement/value now include the requested validation samples, percentiles, wavelength restriction, grey rescaling and magnitude dependence (`refs/montegriffo_2023_gaia_xp.txt:1042–1071`). The 111-SPSS/about-1% claim matches lines 197–202; malformed quotation removed. |
| obs.gaia-xp-ecs-spd-sed-units | accept | Distinct photon/energy reconstruction equations, corrected locator, wavelength limits and overlap weights match `refs/montegriffo_2023_gaia_xp.txt:526–550`. Primary TeX anchor byte-matches `.tex:706`. Malformed quotation removed; unresolved conversion-factor exclusion retained. |
| obs.notation-gaia-sky | accept | Unsupported radial-velocity `units_source` removed; `cm s^-1` retained as project convention. Definitions and epoch match `refs/lindegren_2021_gaia_astrometry.txt:432–437,470–475`; both anchors byte-match. |
| obs.notation-gaia-photometry | accept | Canonical mean-flux symbols corrected, source spellings retained, “unfiltered” removed, and unsupported pupil/response source-unit fields deleted. Definitions match `refs/riello_2021_gaia_photometry.txt:934–940,999–1000,1108–1140,1147–1156`; both anchors byte-match. |

- Mechanical check, run first, exit 0: `493 entries: {'proposed': 52, 'blessed': 350, 'rejected': 91}; 0 errors, 0 warnings`.
- No remaining round-2 notation corrections: source spellings and canonical symbols are distinguished; the per-cm/per-Å distinction is documented.
- XP acceptance covers the restricted statement; the unresolved conversion factor remains excluded from certification.