# Notation additions, lane sky (astrometric symbols owned by obs.notation-gaia-sky, photometric symbols by obs.notation-gaia-photometry)

Dictionary row format: symbol | meaning | units | source spellings (entries) | adopted | note.

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `\alpha_{RA}` | right ascension (ICRS, Gaia-CRF3) of a source at the reference epoch J2016.0 | rad (catalogue: deg) | `alpha` (obs.notation-gaia-sky) | obs.notation-gaia-sky | Bare alpha is the emission angle. |
| `\delta_{Dec}` | declination (ICRS, Gaia-CRF3) of a source at J2016.0 | rad (catalogue: deg) | `delta` (obs.notation-gaia-sky) | obs.notation-gaia-sky | Bare delta is an alias of the emission angle alpha (Pech83). |
| `\varpi_{plx}` | trigonometric parallax at J2016.0 (Gaia astrometric parameter) | rad (catalogue: mas) | `varpi` (obs.notation-gaia-sky) | obs.notation-gaia-sky | Lead suggested bare varpi; subscripted to keep clear of the albedo rows `\varpi_\nu`, `\varpi_\nu^j`. |
| `\mu_{\alpha*}` | proper motion in right ascension times cos delta (pmra) | rad s^-1 (catalogue: mas yr^-1) | `mu_{alpha*} = mu_alpha cos(delta)` (obs.notation-gaia-sky) | obs.notation-gaia-sky | Not the ray cosine mu. |
| `\mu_\delta` | proper motion in declination (pmdec) | rad s^-1 (catalogue: mas yr^-1) | `mu_delta` (obs.notation-gaia-sky) | obs.notation-gaia-sky | Not the ray cosine mu. |
| `v_r` | radial velocity relative to the solar system barycentre | cm s^-1 (catalogue: km s^-1) | `v_r` (obs.notation-gaia-sky) | obs.notation-gaia-sky | |
| `A_{u}` | astronomical unit (in mu_r = v_r varpi / A_u) | cm | `A_u` (obs.notation-gaia-sky) | obs.notation-gaia-sky | Check against const_iau.toml for the project's au symbol before blessing. |
| `G_{Gaia}` | Gaia G-band magnitude | mag | `G` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | Bare G is Newton's constant. |
| `G_{BP}` | Gaia blue-photometer band magnitude | mag | `G_BP` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | |
| `G_{RP}` | Gaia red-photometer band magnitude | mag | `G_RP` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | |
| `n_{p,Gaia}` | mean source flux in a Gaia band in e- s^-1 (not normalised by telescope area); in the XP paper n_p(lambda) is the photon flux distribution | s^-1 (e- s^-1) | `n_p` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | Bare n_p clashes with opacity.notation-particles-plasma (renamed to avoid the check.py conflict). |
| `S_{Gaia}` | Gaia system overall response function (passband) incl. photon-to-electron scaling | e- photon^-1 | `S(lambda)` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | Bare S is the source function. |
| `P_{pupil}` | telescope pupil area in the Gaia flux integral | cm^2 | `P` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | |
| `ZP_{VEG}` | Gaia passband zero point, VEGAMAG, per band | mag | `ZP_{VEG}` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | |
| `ZP_{AB}` | Gaia passband zero point, AB, per band | mag | `ZP_{AB}` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | |
| `\langle F_\lambda\rangle` | synthetic mean energy flux per wavelength unit through a passband | erg s^-1 cm^-2 cm^-1 (printed: W m^-2 nm^-1) | `<f_lambda>` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | Uppercase F because canonical f_nu is the Eddington factor; source spelling lowercase. |
| `\langle F_\nu\rangle` | synthetic mean energy flux per frequency unit through a passband | erg s^-1 cm^-2 Hz^-1 (printed: W m^-2 Hz^-1) | `<f_nu>` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | Uppercase F because canonical f_nu is the Eddington factor; source spelling lowercase. |

Round-1 review (sky-r1) split the original bundled entry: Lindegren 2021 supports only the astrometric rows, so the photometric rows moved to `obs.notation-gaia-photometry` (Riello 2021). Bare aliases `alpha`, `delta`, `varpi`, `G`, `n_p` are no longer registered; the source spelling is kept in `sym_source`. Round 2 (sky-r2) adopted `\langle F_\lambda\rangle`, `\langle F_\nu\rangle` for the synthetic fluxes because `f_nu` is the Eddington factor (collision); source spellings stay lowercase.
