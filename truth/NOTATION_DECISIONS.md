# Notation decisions (canonical)

Project-wide notation for the neutron-star pipeline, decided 2026-10-23 in the notation consolidation pass. **This file is canonical.** `truth/NOTATION.md` is the per-entry rendering produced by `check.py`. This file records the decisions behind it. Every `sym` string in a `kind = "notation"` entry matches exactly one dictionary row below. When you add notation, use the spelling given here, put the paper's spelling in `sym_source`, and add a row here when the quantity is new.

## Conventions

- **Units.** Gaussian CGS throughout (cm, g, s, K, erg, dyn, G, esu). Angles are in rad in code; sources that quote degrees (Ho07 tables) are converted. Photon energies may be quoted in keV: 1 keV = 1.602176634e-9 erg, E = h nu, and I_E = I_nu / h.
- **CODATA.** The `const.*` entries record CODATA 2018 values **in SI, as printed**. Code converts them: kg -> g x 1e3; m -> cm x 1e2; J -> erg x 1e7; C -> esu x 2.99792458e9 (= 10 c[m/s]); J T^-1 -> erg G^-1 x 1e3; W m^-2 K^-4 -> erg cm^-2 s^-1 K^-4 x 1e3; W m^2 -> erg cm^2 s^-1 x 1e11; m K -> cm K x 1e2. `\varepsilon_0` and `\mu_0` have no Gaussian role and are kept only for cross-checks. Formulae with e^2 are Gaussian (no 4 pi epsilon_0). The esu conversion e[esu] = 10 c[m/s] e[C] is the project's definition of the statcoulomb (no local source prints it exactly; see const.esu-from-coulomb, rejected r2, work order: Jackson units appendix). With it, e[esu] = 4.803204712570e-10 while sqrt(alpha hbar c) = 4.803204713885e-10 with the exact hbar = h/(2 pi) and CODATA 2018 alpha; the relative gap is 2.74e-10 in e (5.5e-10 in e^2), which equals half of mu_0/(4 pi 1e-7) - 1 = 5.44e-10: in the 2019 SI mu_0 is no longer exactly 4 pi 1e-7, so the Gaussian identity e^2 = alpha hbar c holds only to that level. An earlier claim of a 3.3e-11 gap used the truncated printed hbar and was wrong; corrected 2026-10-23 (40-digit recomputation). Code uses the conversion for e and never asserts e^2 = alpha hbar c beyond 1e-9. Decision (pipeline-r2, 2026-10-23): code computes the Gaussian e as sqrt(alpha_f hbar c) from CGS alpha_f, hbar, c (the self-consistent Gaussian definition); the SI bridge 10 c e_SI is a cross-check agreeing to 2.7e-10 and is not a runtime input. Spectral Jacobians: I_E = I_nu/h per erg, multiply by 1.602176634e-9 erg/keV for per-keV units; I_lambda = I_nu c/lambda^2 with 1 nm = 1e-7 cm and 1 A = 1e-8 cm.
- **Masses.** `m_u` is the CODATA atomic mass constant. `m_u^{BSk}` = 1.66e-24 g is used only inside the BSk fits. `m_H` = m_p + m_e is PC03's definition; McPHAC's "hydrogen mass" is not silently identified with it.
- **Logarithms.** Write `ln` for the natural log and `log10` explicitly. Only the BSk entries use log = log10.
- **Observer frame.** Superscript infinity means "measured at infinity": `T_{eff}^\infty`, `L_{ph}^\infty`, `R^\infty`, `\nu^\infty`, `I^\infty`. Redshift is `z_g` with 1 + z_g = (1 - r_g/R)^(-1/2). Subscript infinity is never used, which avoids clashing with the Rydberg constant (`R_{Ry}`). Camera-frame quantities carry the subscript `cam`.
- **Directions.** `\theta` is the angle between the ray and the outward surface normal, and `\mu` = cos theta. McPHAC and VAL06 use a half-range mu > 0 with explicit +-mu for outgoing and ingoing rays; SPW09 uses a signed mu. The project convention is the signed mu = cos theta in [-1, 1]; half-range sources are converted at the interface. `\alpha` is the GR emission angle; for an emergent ray theta = alpha. `\varphi_k` is the photon azimuth about the normal, measured from the (B, n) plane. `\theta_{kB}` is the k-B angle, `\theta_{Bn}` is the B-normal angle, and `\theta_B` is the magnetic colatitude.
- **Indices.** `j` is the normal mode: 1 = X (extraordinary), 2 = O (ordinary); a scattering pair is j -> j'. `\alpha_{pol}` in {-1, 0, +1} is the cyclic polarisation index; in the subscripts of opacity quantities the dummy alpha means alpha_pol. `l` in {1, 2} is the linear basis (in or perpendicular to the plane of incidence) and appears only in the condensed-surface reflectivities. `q` is the Gauss-Legendre node index. `n_L` is the Landau level. Spatial tensor indices are never i, j: the conduction tensor uses `K_{ab}^{cond}` with a, b in {x, z}.
- **Frames for tensors.** Cyclic components `e_{j,\alpha}` and the dielectric elements `\varepsilon_{\perp,\wedge,\parallel}` use the frame z || B. `K_j` and `K_{z,j}` use the frame z || k with B in the x-z plane. Converting between the two needs an explicit rotation.
- **Depth.** Column density is `m` (g cm^-2), with dm = -rho dz. McPHAC and the pipeline entries write it y.
- **Spelling rules.** Symbols are LaTeX, and multi-letter sub- and superscripts are braced (`T_{eff}`). Text labels are plain-letter subscripts (pix, vis, cap, geo, rel). `sym` contains no Unicode Greek; Unicode survives only in `sym_source`. Compound rows (`\mu_q, a_q`) are kept where the source bundles symbols.

## Database fields added by this pass (notation entries only)

- Per symbol: `sym_source` holds the paper's spelling when `sym` was changed. `aliases` lists the other spellings of the same quantity (from source rows and from text). `meaning_source` and `units_source` hold the pre-normalisation text wherever a canonical meaning or unit was imposed. check.py compares meanings literally, so every row that shares a `sym` carries the identical canonical meaning.
- Per entry: `conflicts_with` is the union of the collision groups below that the entry belongs to. `adopted = true` means the entry is the adopted party in at least one group. `adopted_syms` names those groups, and `conflict_notes` repeats the resolutions. check.py reads `adopted` at entry level, so use `adopted_syms` for the answer per symbol. No two entries now share a `sym` with different meanings. If a future literal clash appears, check.py will again require exactly one adopted party.

## Retired bare symbols

These letters had too many meanings, so they are never used without decoration: `u` (now `m_u`, `x_g`, `u_{geo}`, `u_{kT}`, `u_e`, `u_i`, `u_r`, `u^\mu_{obs}`), `L` (now `L_{ph}`, `L_{sc}`), `\xi` (now `\xi_\rho`, `\xi_{sky}`, `\xi_{vac}`), `n` (now `n_L`, `n_{fm}`, `n_e`, `n_b`), `H` (now `\mathcal{H}`, `H_\nu`, `H_0`, `H_\rho`) and `S` (now `S_\nu`, `S_{vis}`, `S_{stim}`). Bare `y` and `z` are CIE chromaticity only, and bare `Z` is the CIE tristimulus value only.

## Collisions resolved (one source spelling, several meanings)

| source spelling | adopted meaning (entry) | parties | resolution |
|---|---|---|---|
| `alpha` | gr.notation-geometry | gr.notation-pechenick, const.notation-atomic, opacity.notation-particles-plasma, opacity.notation-modes-polarization, obs.notation-rxj1856-photometry, rt.notation-spw09, surface.notation-geometry, rt.notation-cond-energies-geometry | Bel02 emission angle keeps \alpha. Renamed: Pech83 polar-cap radius -> \theta_{cap}; CODATA/PLCH04 fine-structure -> \alpha_f; PC03/PLCH04 cyclic polarisation index -> \alpha_{pol}; Kaplan11 optical index -> \alpha_{opt}; SPW09 B-k angle -> \theta_{kB}; SPW09 albedo -> \varpi_\nu^j; GH83 pole-to-sight angle -> \psi_{pole}; Pot12 incident/reflected k-B angles -> \theta_{kB}^{(i)}, \theta_{kB}^{(r)}, \theta_{kB}^{(min)}. |
| `lambda` | colour.notation-cie-cmf | gr.notation-geometry, gr.notation-pechenick | Wavelength keeps \lambda. Bel02/Pech83 affine parameter -> \lambda_{aff}. |
| `phi` | gr.notation-geometry | gr.notation-pechenick, rt.notation-spw09, rt.notation-cond-energies-geometry, surface.notation-geometry | Bel02 surface-element azimuth about the line of sight keeps \phi. Pech83 orbit-plane azimuth -> \phi_{orb}; SPW09/Pot12 photon azimuth about the normal -> \varphi_k; GH83 rotation phase -> \phi_{rot}. |
| `theta` | rt.notation-spw09 | gr.notation-geometry, gr.notation-pechenick, surface.notation-geometry, rt.notation-cond-energies-geometry | RT ray angle to the surface normal keeps \theta (mu = cos theta). Schwarzschild polar coordinate (Bel02, Pech83) -> \theta_{Sch}; GH83 line-of-sight colatitude -> \theta_{los}; Pot12 theta_k merged into \theta. |
| `theta_B` | surface.notation-geometry | opacity.notation-modes-polarization, rt.notation-cond-energies-geometry, rt.notation-spw09, obs.notation-rxj1856-spectral | Magnetic colatitude keeps \theta_B. PC03/PLCH04 k-B angle -> \theta_{kB}; SPW09/Pot12/Ho07 field-normal angle -> \theta_{Bn}. |
| `mu` | rt.notation-mcphac | gr.notation-geometry, eos.notation-electron-gas, opacity.notation-particles-plasma | RT direction cosine keeps \mu. Bel02 mu = cos psi -> \cos\psi; PY01 chemical potential -> \mu_e; PC03 reduced mass -> \mu_{red}. |
| `b` | gr.notation-geometry | eos.notation-electron-gas, opacity.notation-fields | Impact parameter keeps b. Relativistic field parameter B/B_Q -> b_{rel}. |
| `u` | none: bare symbol retired | const.notation-particles, gr.notation-geometry, gr.notation-pechenick, opacity.notation-photon-cyclotron | Bare u retired. CODATA atomic mass unit -> m_u; Bel02 compactness -> x_g; Pech83 M/r -> u_{geo}; PC03 hbar omega/kT -> u_{kT}. |
| `R_inf` | surface.notation-observer-frame | const.notation-atomic, gr.notation-potekhin, eos.notation-structure, obs.notation-rxj1856-spectral | Apparent stellar radius is R^\infty (one spelling). CODATA Rydberg R_inf -> R_{Ry}; Ho07 emission-region apparent radius -> R_{em}^\infty. |
| `h` | const.notation-fundamental | eos.notation-bsk-fit | Planck constant keeps h. BSk enthalpy per baryon -> h_b (surface value h_{b,s}). |
| `B` | opacity.notation-fields | gr.notation-pechenick | Magnetic field keeps B. Pech83 effective potential B^-2(r) -> \mathcal{B}^{-2}(r). |
| `rho` | rt.notation-spw09 | rt.notation-mcphac | Mass density keeps \rho. McPHAC albedo rho_nu -> \varpi_\nu. |
| `k` | rt.notation-mcphac | rt.notation-spw09, const.notation-fundamental | McPHAC k_nu (B = 0 total opacity) keeps k_\nu. SPW09 k_nu^i (mode absorption) -> \kappa^a_j; CODATA k -> k_B. |
| `sigma` | opacity.notation-opacity-symbols | const.notation-fundamental, rt.notation-spw09 | sigma is a cross section (cm^2). CODATA Stefan-Boltzmann sigma -> \sigma_{SB}; SPW09 scattering opacity sigma_nu^i -> \kappa^s_j, differential sigma_nu^ij -> d\kappa^s_{ji}/d\Omega. |
| `tau` | rt.notation-mcphac | rt.notation-val | Monochromatic optical depth keeps \tau_\nu. VAL06 Thomson depth tau -> \tau_T. |
| `j` | rt.notation-spw09 | rt.notation-mcphac, rt.notation-cond-emissivity-fit | j is the normal-mode index (1 = X, 2 = O). McPHAC angle-node index j -> q; Pot12 linear-polarisation index j -> l (m -> l'), Pot12 mode index M -> j. |
| `M` | gr.notation-geometry | gr.notation-pechenick, rt.notation-cond-emissivity-fit | M is the gravitational mass in g. Pech83 geometric-unit M -> GM/c^2; Pot12 mode label M -> j. |
| `Z` | colour.notation-cie-xyz | opacity.notation-particles-plasma, rt.notation-cond-matter-dielectric | Bare Z is the CIE tristimulus value. Ion charge number -> Z_{ion} (mass number A -> A_{ion}). |
| `L` | none: bare symbol retired | eos.notation-structure, colour.notation-reinhard-luminance | Bare L retired. PPP15 local luminosity -> L_{ph}; Reinhard scaled luminance -> L_{sc}(x_{pix},y_{pix}). |
| `y` | colour.notation-chromaticity | colour.notation-reinhard-parameters, rt.notation-mcphac | Chromaticity keeps x, y, z. Reinhard pixel position -> (x_{pix},y_{pix}); McPHAC column depth y -> m. |
| `xi` | none: bare symbol retired | eos.notation-bsk-fit, gr.notation-pechenick, opacity.notation-vacuum | Bare xi retired. BSk log10 rho -> \xi_\rho; Pech83 sky angle -> \xi_{sky}; PLCH04 1/b -> \xi_{vac}. |
| `eta` | rt.notation-spw09 | eos.notation-bsk-fit, rt.notation-cond-matter-dielectric, opacity.notation-modes-polarization | eta is an emissivity (\eta_\nu^j). BSk exp(H)-1 -> \eta_{\mathcal{H}}; vAL05/PLCH04 longitudinal dielectric element -> \varepsilon_\parallel. |
| `g` | rt.notation-spw09 | opacity.notation-modes-polarization, rt.notation-cond-matter-dielectric | g is surface gravity. vAL05/PLCH04 gyrotropic dielectric element g -> \varepsilon_\wedge (PHC16 spelling); transverse element epsilon -> \varepsilon_\perp. |
| `c_i` | const.notation-fundamental | eos.notation-bsk-fit | c_1, c_2 are the radiation constants. BSk fit coefficients a_i, b_i, c_i -> a_i^{BSk}, b_i^{BSk}, c_i^{BSk}. |
| `n` | none: bare symbol retired | eos.notation-bsk-fit, eos.notation-electron-gas | Bare n retired. BSk n_b/fm^-3 -> n_{fm}; Landau level index -> n_L. |
| `S` | rt.notation-mcphac | gr.notation-geometry, rt.notation-vacuum | S_nu is the source function. Bel02 visible area S_v -> S_{vis}; VAL06 stimulated-emission factor S -> S_{stim}. |
| `nu` | gr.notation-radiation | obs.notation-rxj1856-timing | nu is a photon frequency. Spin frequency -> \nu_{spin}. |
| `beta` | gr.notation-geometry | opacity.notation-photon-cyclotron, opacity.notation-modes-polarization, surface.notation-geometry | Bel02 bending angle keeps \beta. PC03 beta_e, beta_p -> \zeta_e, \zeta_p (PHC16 spelling); PLCH04 mode parameter -> \beta_{pol}; GH83 beta (magnetic colatitude) -> \theta_B. |
| `H` | rt.notation-spw09 | eos.notation-bsk-fit | H_nu is the Eddington flux. BSk pseudo-enthalpy H -> \mathcal{H}. |

## Dictionary

Columns: project symbol (the `sym` string in the database, LaTeX); meaning; units (CGS unless stated); source spellings with the notation entries that use them; adopted notation entry; note. Rows are generated from the database, one per distinct `sym`.

### Constants (CODATA 2018; SI as printed, CGS in code)

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `h` | Planck constant | erg s | `h` (const.notation-fundamental) | const.notation-fundamental | SI in CODATA source (J Hz^-1). Bare h is reserved for Planck; BSk enthalpy per baryon renamed h_b. (r1 fix: CODATA spells the symbol h, not h_P; SI value h = 6.626 070 15e-34 J Hz^-1, Table XXXI.) |
| `\hbar` | reduced Planck constant h/(2 pi) | erg s | `ℏ` (const.notation-fundamental); `hbar` (opacity.notation-particles-plasma) | const.notation-fundamental |  |
| `e` | elementary charge | esu (statC) | `e` (const.notation-fundamental, opacity.notation-particles-plasma) | opacity.notation-particles-plasma | Gaussian CGS everywhere in code; CODATA prints C: e[esu] = e[C] x 10 c[m s^-1]. Not the polarisation vector (bold e_j). (r1 reject upheld 2026-10-23: the exact conversion e[esu] = 10 c[m/s] e[C] disagrees with sqrt(alpha hbar c) by 2.74e-10 in e once the exact hbar = h/(2 pi) is used; the earlier 3.3e-11 figure used the truncated printed hbar. The conversion is kept as the project's definition; e^2 = alpha hbar c holds only to 5.5e-10 in e^2 because mu_0 is not exactly 4 pi 1e-7 in the 2019 SI.) Decision (pipeline-r2, 2026-10-23): code computes the Gaussian e as sqrt(alpha_f hbar c) from CGS alpha_f, hbar, c (the self-consistent Gaussian definition); the SI bridge 10 c e_SI is a cross-check agreeing to 2.7e-10 and is not a runtime input. |
| `k_B` | Boltzmann constant | erg K^-1 | `k` (const.notation-fundamental); `k_B` (opacity.notation-particles-plasma) | opacity.notation-particles-plasma |  |
| `N_A` | Avogadro constant | mol^-1 | `N_A` (const.notation-fundamental) | const.notation-fundamental | (r1 fix applied: locator widened to 'Table XXXI, pp. 025010-45, 025010-49' in const.notation-fundamental (truth/entries/const.toml:640); N_A at refs/tiesinga_2021_codata2018.txt L3730, exact value, p. 49; per-symbol source field absent in check.py, so recorded in entry notes.) |
| `\sigma_{SB}` | Stefan-Boltzmann constant | erg cm^-2 s^-1 K^-4 | `σ` (const.notation-fundamental); `σ_sb` (rt.notation-val); `σ_SB` (surface.notation-observer-frame); text/other: `sigma_SB` | surface.notation-observer-frame |  |
| `c_1` | first radiation constant, 2 pi h c^2 | erg cm^2 s^-1 | `c_1` (const.notation-fundamental) | const.notation-fundamental | SI (W m^2) in source. (r1 fix applied: Table XXXI p. 025010-49, tiesinga_2021_codata2018.txt L3765 "first radiation constant 2πhc2"; SI W m^2 before conversion; locator and notes in truth/entries/const.toml:640.) |
| `c_2` | second radiation constant, hc/k | cm K | `c_2` (const.notation-fundamental) | const.notation-fundamental | SI (m K) in source. (r1 fix applied: Table XXXI p. 025010-49, tiesinga_2021_codata2018.txt L3766 "second radiation constant hc=k"; locator and notes in truth/entries/const.toml:640.) |
| `\varepsilon_0` | vacuum electric permittivity | F m^-1 (SI only) | `ε_0` (const.notation-fundamental) | const.notation-fundamental | Not used in Gaussian CGS formulae; kept for CODATA cross-checks. |
| `\mu_0` | vacuum magnetic permeability | N A^-2 (SI only) | `μ_0` (const.notation-fundamental) | const.notation-fundamental | Not used in Gaussian CGS formulae. Not the ray cosine mu. |
| `\mathrm{eV}` | electron volt, an energy unit: 1 eV = 1.602176634e-19 J | erg | `eV` (const.notation-fundamental) | const.notation-fundamental | 1 eV = 1.602176634e-12 erg; 1 keV = 1.602176634e-9 erg. |
| `m_e` | electron mass | g | `m_e` (const.notation-particles, opacity.notation-particles-plasma) | const.notation-particles | CODATA kg x 1e3. |
| `m_p` | proton mass | g | `m_p` (const.notation-particles, opacity.notation-particles-plasma) | const.notation-particles | CODATA kg x 1e3. |
| `m_n` | neutron mass | g | `m_n` (const.notation-particles) | const.notation-particles |  |
| `m_u` | atomic mass constant m_u (= unified atomic mass unit u) | g | `u` (const.notation-particles) | const.notation-particles | CODATA atomic mass constant. Never substitute the BSk rounded m_u^{BSk} = 1.66e-24 g outside the BSk fits. |
| `m_p/m_e` | proton-electron mass ratio | 1 (SI as printed) | `m_p/m_e` (const.notation-particles) | const.notation-particles |  |
| `\alpha_f` | fine-structure constant e^2/(hbar c) | 1 (SI as printed) | `α` (const.notation-atomic); `alpha_f` (opacity.notation-particles-plasma); text/other: `alpha` | const.notation-atomic | (r1 fix applied: Gaussian definition printed in potekhin_lai_chabrier_ho_2004_polarization L544 "where αf = e2 /~c = 1/137.036 is the fine-structure constant"; CODATA SI formula at L3481; recorded in const.notation-atomic notes, truth/entries/const.toml:696; no per-symbol source field in check.py.) |
| `r_e` | classical electron radius, alpha^2 a_0 | cm | `r_e` (const.notation-atomic) | const.notation-atomic |  |
| `\sigma_T` | Thomson cross section (8 pi/3) r_e^2, r_e = e^2/(m_e c^2) | cm^2 | `σ_e` (const.notation-atomic); `sigma_T` (opacity.notation-opacity-symbols) | const.notation-atomic |  |
| `a_0` | Bohr radius | cm | `a_0` (const.notation-atomic) | const.notation-atomic |  |
| `R_{Ry}` | Rydberg constant | cm^-1 | `R_∞` (const.notation-atomic) | const.notation-atomic | CODATA writes R_∞; renamed so that the superscript/subscript infinity is never the Rydberg constant. |
| `cR_{Ry}` | Rydberg frequency | Hz | `cR_∞` (const.notation-atomic) | const.notation-atomic |  |
| `\mathrm{Ry}` | Rydberg energy | erg | `hcR_∞` (const.notation-atomic) | const.notation-atomic | Rydberg energy hcR_∞ (McPHAC gamma^2 = Z^2 Ry/(k_B T)). |
| `E_h` | Hartree energy, 2 hcR_inf | erg | `E_h` (const.notation-atomic) | const.notation-atomic |  |
| `\mu_B` | Bohr magneton: SI e hbar/(2 m_e); Gaussian e hbar/(2 m_e c) = E_ce/(2B) | erg G^-1 | `μ_B` (const.notation-atomic) | const.notation-atomic | SI J T^-1 in source; x 1e3. (r1 fix: CODATA prints the SI form e hbar/(2 m_e); Gaussian-unit use needs the extra 1/c, equivalently E_ce/(2B).) |
| `u_r` | relative standard uncertainty | 1 (SI as printed) | `u_r` (const.notation-atomic) | const.notation-atomic |  |
| `G` | gravitational constant | cm^3 g^-1 s^-2 | `G` (eos.notation-structure) | eos.notation-structure |  |
| `c` | speed of light | cm s^-1 | `c` (eos.notation-structure, opacity.notation-particles-plasma) | eos.notation-structure |  |
| `m_H` | hydrogen atom mass m_p + m_e | g | `m_H` (opacity.notation-particles-plasma) | opacity.notation-particles-plasma |  |
| `\mathrm{au}` | astronomical unit, exactly 149 597 870 700 m (IAU 2012 B2) | cm | `au` (const.au); `A_u` (obs.notation-gaia-sky, obs.gaia-perspective-acceleration, obs.agis-astrometric-parameter-set, obs.agis-standard-model-direction) | const.au | (M1 2026-10-23) SI (m) as printed; x 1e2 in code. A unit symbol like `\mathrm{eV}`, `\mathrm{pc}`. Replaces the sky lane's `A_{u}` (entry edit pending, truth/harvest/notation_merge_report.md). |
| `\mathrm{pc}` | parsec, exactly 648000/pi au (IAU 2015 B2) | cm | `pc` (const.parsec) | const.parsec | (M1 2026-10-23) SI (m) as printed; defined through `\mathrm{au}`. Obs tables quoting pc convert through this row (answers the r1 cross-row finding on the parsec). |
| `(GM)^{N}_\odot` | nominal solar mass parameter, 1.3271244e20 m^3 s^-2 (conversion constant, IAU 2015 B3) | cm^3 s^-2 | `\mathcal{(GM)}^{\mathrm N}_\odot` (const.gm-sun-nominal, const.gm-sun-nominal-cgs, const.m-sun-nominal-cgs) | const.gm-sun-nominal | (M1 2026-10-23) SI as printed; const.gm-sun-nominal-cgs gives the cm^3 s^-2 value. A solar-mass unit follows only through G (const.m-sun-nominal-cgs); no row for it yet. |

### GR exterior, ray geometry and stellar structure

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `r` | radial coordinate measured from the stellar centre (areal radius) | cm | `r` (eos.notation-structure, gr.notation-geometry, gr.notation-pechenick, gr.notation-potekhin) | gr.notation-geometry |  |
| `R` | stellar (circumferential) radius | cm | `R` (eos.notation-structure, gr.notation-geometry, gr.notation-pechenick, gr.notation-potekhin); `R_NS` (rt.notation-spw09) | gr.notation-geometry |  |
| `M_r` | gravitational mass inside a sphere of radius r, M_r = 4 pi int_0^r rho r'^2 dr' | g | `M_r` (eos.notation-structure, gr.notation-potekhin) | gr.notation-potekhin |  |
| `M` | total gravitational stellar mass, M = M_R | g | `M` (eos.notation-structure, gr.notation-geometry, gr.notation-potekhin); `M_NS` (rt.notation-spw09) | gr.notation-geometry | Mass in g. Pechenick geometric-unit M is GM/c^2. |
| `r_g` | Schwarzschild radius 2GM/c^2 of the star | cm | `r_g` (eos.notation-structure, gr.notation-geometry, gr.notation-potekhin); `R_S` (rt.notation-spw09) | gr.notation-geometry |  |
| `\Phi` | metric function Phi(r): g_tt = -e^(2 Phi) c^2 | 1 | `Phi(r)` (eos.notation-structure); `Phi` (gr.notation-potekhin) | gr.notation-potekhin |  |
| `\Lambda` | metric function Lambda(r): g_rr = e^(2 Lambda) | 1 | `Lambda(r)` (eos.notation-structure); `Lambda` (gr.notation-potekhin) | gr.notation-potekhin | Not the Lambda operator or Coulomb logarithm (always subscripted). The project-derived Schwarzschild identity e^(-2 Lambda) = 1 - 2GM_r/(c^2 r) is not asserted as a source statement: PPP15 eq. 1 prints e^(2 Lambda) = 1 - 2GM_r/(c^2 r), which is inconsistent with g_rr = e^(2 Lambda) (r1 fix applied; see gr r1 review). |
| `P` | pressure | dyn cm^-2 | `P` (eos.notation-structure) | eos.notation-structure |  |
| `\rho` | mass density (TOV/EOS: gravitational mass density = energy density/c^2; atmosphere/opacity: plasma mass density, equal to it up to O(k_B T/(m_u c^2))) | g cm^-3 | `rho` (eos.notation-structure); `rho_gas` (opacity.notation-particles-plasma); `ρ` (rt.notation-spw09) | rt.notation-spw09 | (disputed: r1 rejected merging the gravitational mass density (TOV) with the plasma mass density; no supplied source establishes the claimed O(k_B T/(m_u c^2)) equivalence (potekhin_pons_page_2015_thermal:214-215; suleimanov_potekhin_werner_2009_mag_atm:140-143). Row kept as one symbol with the two meanings; use the atmosphere meaning only inside atmosphere/opacity code and the TOV meaning only in structure code.) |
| `g` | surface gravity (local proper) g = GM/(R^2 sqrt(1 - r_g/R)), taken uniform through the atmosphere | cm s^-2 | `g` (eos.notation-structure); `g_s` (rt.notation-mcphac); `g_surf` (rt.notation-spw09) | rt.notation-spw09 |  |
| `t` | Schwarzschild coordinate time | s | `t` (gr.notation-geometry) | gr.notation-geometry |  |
| `\theta_{Sch}` | Schwarzschild polar coordinate (orbit plane taken as theta = 90 deg) | rad | `theta` (gr.notation-geometry, gr.notation-pechenick) | gr.notation-geometry |  |
| `\psi` | polar angle of the photon position in its orbit plane, measured from the escape direction (psi = 0 at infinity); for a surface point, the polar angle of the emission point measured from the direction to the observer | rad | `psi` (gr.notation-geometry); `Delta phi` (gr.notation-pechenick); text/other: `Δφ` | gr.notation-geometry | (r1 fix: Pech83 Delta phi_orb for a distant observer is psi evaluated at emission; (Delta phi)_max is psi_max = psi(b_max).) |
| `\alpha` | emission angle: angle between the photon direction and the local radial direction (surface normal), measured by the local static observer at the emission point | rad | `alpha` (gr.notation-geometry); `delta` (gr.notation-pechenick) | gr.notation-geometry | For an emergent ray theta (RT) = alpha and mu = cos alpha. Bare alpha means only this. (r1 fix blocked: the alpha/delta alias merge changes the meaning of theta across rt.notation-spw09 and gr.notation-geometry, and gr.notation-geometry is blessed; needs an owner decision, see Remaining doubts. No edit.) |
| `\beta` | bending angle beta = psi - alpha | rad | `beta` (gr.notation-geometry) | gr.notation-geometry | Bending angle only; beta_pol, zeta_e etc. are distinct. |
| `b` | impact parameter of the photon orbit | cm | `b (impact)` (gr.notation-geometry, gr.notation-pechenick) | gr.notation-geometry | Impact parameter only; the relativistic field parameter is b_rel. |
| `\cos\psi` | mu = cos psi | 1 | `mu (cos psi)` (gr.notation-geometry) | gr.notation-geometry | Bel02 writes mu; renamed so that bare mu is the RT direction cosine. |
| `\cos\psi_{vis}` | critical value of mu = cos psi at which the emission is tangential (alpha = 90 deg); surface points with mu < mu_v are invisible (single-image regime R > 2 r_g only; compact stars with R < 2 r_g or R < 3 GM/c^2 show multiple images) | 1 | `mu_v` (gr.notation-geometry) | gr.notation-geometry |  |
| `x_g` | compactness r_g/R = 2GM/(c^2 R) at the stellar surface | 1 | `u (r_g/R)` (gr.notation-geometry); `x_g` (surface.notation-observer-frame); text/other: `u` | surface.notation-observer-frame |  |
| `\lambda_{aff}` | affine parameter along the photon orbit (normalisation source-specific: Bel02 lambda -> r as r -> infinity; Pech83 eq. 2.2) | cm | `lambda` (gr.notation-geometry, gr.notation-pechenick) | gr.notation-geometry | (disputed: r1 rejected a canonical affine-parameter normalisation; each source keeps its own (Bel02 lambda -> r as r -> infinity; Pech83 eq. 2.2). Do not convert between them until one coordinate convention is adopted.) |
| `k^\mu` | photon tangent four-vector dx^mu/dlambda_aff (normalisation source-specific: Bel02 u_t = 1; Pech83 eq. 2.2) | s cm^-1 for u^t; 1 for u^r; cm^-1 for u^psi (t in s, lambda in cm; source normalisation u_t = 1 is source-specific, conversion to the adopted metric convention unspecified) | `u^k` (gr.notation-geometry); `v^mu` (gr.notation-pechenick) | gr.notation-geometry | (disputed: r1 rejected assigning CGS component units while the conversion between the Bel02 (u_t = 1) and Pech83 (eq. 2.2) normalisations is unspecified. The two source tangents dx/dlambda are separate; their common project normalisation is not established.) |
| `D` | distance from the observer to the star (D >> R) | cm | `D` (gr.notation-geometry); `d` (obs.notation-rxj1856-spectral) | gr.notation-geometry | Obs tables quote pc (1 pc = 3.0857e18 cm). (r1 fix applied: the uncited parsec conversion was removed from this note; a finite camera radius is a separate symbol.) |
| `\phi` | azimuthal angle of a surface element about the observer direction D | rad | `phi` (gr.notation-geometry) | gr.notation-geometry | Azimuth of a surface element about the line of sight only. |
| `S_{vis}` | visible area of the stellar surface | cm^2 | `S_v` (gr.notation-geometry) | gr.notation-geometry |  |
| `I_0` | surface radiation intensity at the emission point as a function of emission angle alpha (emitter frame); the paper relates it to the observed I through I/nu^4 = I_0/nu_0^4 | erg s^-1 cm^-2 sr^-1 (inferred frequency-integrated, from the I/nu^4 scaling) | `I_0` (gr.notation-radiation) | gr.notation-radiation |  |
| `I^\infty` | observed intensity at the observer (at infinity) corresponding to I_0 | erg s^-1 cm^-2 sr^-1 (inferred, as I_0) | `I` (gr.notation-radiation) | gr.notation-radiation |  |
| `\nu^\infty` | photon frequency measured at the observer (at infinity) | Hz | `nu (observer)` (gr.notation-radiation) | gr.notation-radiation |  |
| `\nu` | photon frequency in the local rest frame of the emitting plasma (static frame at the surface for a non-rotating star) | Hz | `nu_0` (gr.notation-radiation); `nu (I_nu)` (gr.notation-pechenick); `nu_loc` (opacity.notation-photon-cyclotron); text/other: `nu` | gr.notation-radiation | Observer-frame frequency is nu^inf (camera frame nu_cam). Spin frequency is nu_spin. Collision frequencies carry sub/superscripts. (disputed: r1 rejected forcing the frame-generic Pech83 nu and the distant Bel02 nu into the emitting-plasma frame. Survives: Bel02 static local nu_0, the distant nu^inf, and the Pech83 frame-dependent nu, to be labelled separately before any rotating-star work.) |
| `dS` | proper surface element of the stellar surface, dS = R^2 d(cos psi) d phi | cm^2 | `dS` (gr.notation-radiation); `dΣ` (surface.notation-envelope-fit) | gr.notation-radiation |  |
| `d\Omega` | solid angle subtended by dS on the observer sky, dOmega = b db dphi / D^2 | sr | `dOmega` (gr.notation-radiation) | gr.notation-radiation |  |
| `dF` | observed flux from the surface element dS, dF = I dOmega | erg s^-1 cm^-2 (frequency-integrated) | `dF` (gr.notation-radiation) | gr.notation-radiation |  |
| `GM/c^2` | gravitating mass M of the star in geometric units G = c = 1, i.e. the length GM/c^2 | cm | `M (G=c=1)` (gr.notation-pechenick) | gr.notation-pechenick | Pechenick geometric-unit mass (= r_g/2); convert on use. |
| `r_0` | radius of the stationary observer, r_0 -> infinity | cm | `r_0` (gr.notation-pechenick) | gr.notation-pechenick |  |
| `b_{max}` | largest impact parameter of a photon that can leave the surface r = R, b_max = R (1 - r_g/R)^(-1/2) (project units; Pech83 writes R(1-2M/R)^(-1/2) with geometric M); valid for R > 3GM/c^2 | cm | `b_max` (gr.notation-pechenick) | gr.notation-pechenick |  |
| `\phi_{orb}` | azimuthal angle of the photon in its orbit plane (eq. 2.2) | rad | `phi` (gr.notation-pechenick) | gr.notation-pechenick |  |
| `u_{geo}` | dimensionless inverse radius u_geo = GM/(c^2 r) (Pech83 u = M/r with geometric M; running variable of the orbit equation) | 1 | `u (M/r)` (gr.notation-pechenick) | gr.notation-pechenick |  |
| `u_{b,geo}` | u_b,geo = GM/(c^2 b) (Pech83 u_b = M/b with geometric M) | 1 | `u_b` (gr.notation-pechenick) | gr.notation-pechenick |  |
| `\xi_{sky}` | angle in the observer sky between a photon and the radial photon direction (the static observer sees the radial photon along the line to the star centre) | rad | `xi (sky)` (gr.notation-pechenick); text/other: `theta'`, `theta_prime` | gr.notation-pechenick | (r1 fix: Pech83 theta' is an alias of this angle; units rad.) |
| `\psi_{max}` | value of psi at emission for b = b_max (Pech83 (Delta phi)_max) | rad | `(Delta phi)_max` (gr.notation-pechenick) | gr.notation-pechenick |  |
| `u^\mu_{obs}` | four-velocity of the stationary (static) observer | u^t dimensionless (t in s), u^r cm s^-1, u^theta and u^phi s^-1 (coordinate components; static observer: u^r = u^theta = u^phi = 0) | `u^mu` (gr.notation-pechenick) | gr.notation-pechenick |  |
| `\mathcal{B}^{-2}(r)` | B^(-2)(r) = (1 - r_g/r)/r^2 (Pech83 (1-2M/r)/r^2 with geometric M), effective potential of the photon orbit | cm^-2 | `B^-2(r)` (gr.notation-pechenick) | gr.notation-pechenick | Not a magnetic field. |
| `\theta_{cap}` | polar-cap angular radius in Sections II-III (Pech83 alpha for the cap; distinct from the Appendix photon-cone half-angle alpha) | rad | `alpha` (gr.notation-pechenick) | gr.notation-pechenick | Pech83 Sects. II-III; the Appendix reuses alpha for a photon-cone half-angle. |
| `I_\nu` | intensity per unit frequency interval | erg s^-1 cm^-2 Hz^-1 sr^-1 | `I_nu` (gr.notation-pechenick) | gr.notation-pechenick |  |
| `\tilde E` | specific energy (per unit rest mass) of a timelike geodesic, conserved: -u_t in geometric units, (1 - r_g/r) dt/dtau in CGS | 1 | `E` (gr.timelike-constants-schwarzschild, gr.timelike-radial-equation-schwarzschild, gr.timelike-dt-dtau-schwarzschild, gr.timelike-periastron-cubic, gr.timelike-angular-momentum-from-periastron, gr.timelike-unbound-condition, gr.timelike-speed-at-infinity, gr.timelike-proper-time-quadrature; DA09 eqs. 6-10) | gr.timelike-constants-schwarzschild | (M1 2026-10-23) Per unit rest mass, from the norm -delta_1 = -1 of dx/dlambda. DA09 Sect. 3 also writes gamma = E/m; the two agree only for m = 1. Camera value `\tilde E_{cam}`. Unbound (hyperbolic) iff `\tilde E` > 1. Bare E stays the photon energy; `\tilde E_{C}`, `\tilde E_{pe}` are distinct strings. |
| `\tilde L` | specific angular momentum (per unit rest mass) of a timelike geodesic about the orbit-plane normal, conserved: u_phi = r^2 dphi_Sch/dtau_prop | cm^2 s^-1 (geometric: cm) | `L_z` (gr.timelike-constants-schwarzschild, gr.timelike-radial-equation-schwarzschild, gr.timelike-dphi-dtau-schwarzschild, gr.timelike-periastron-cubic, gr.timelike-angular-momentum-from-periastron, gr.timelike-proper-time-quadrature; DA09, equals the total angular momentum for the equatorial orbit, Carter Q = 0) | gr.timelike-constants-schwarzschild | (M1 2026-10-23) Geometric L_z = `\tilde L`/c. Camera value `\tilde L_{cam}`. Sign = orbit orientation, a free choice. Bare L stays retired. |
| `\tau_{prop}` | proper time along a timelike worldline, dtau^2 = -ds^2/c^2 (signature -,+,+,+) | s | `\tau` (gr.timelike-constants-schwarzschild, gr.timelike-dt-dtau-schwarzschild, gr.timelike-dphi-dtau-schwarzschild, gr.timelike-proper-time-quadrature); `\lambda` (DA09 affine parameter, equal to proper time for delta_1 = 1; geometric units, a length c tau) | gr.timelike-constants-schwarzschild | (M1 2026-10-23) Bare tau is taken by optical depths (`\tau_\nu`, `\tau_T`). On the camera worldline this is `\tau_{cam}`. The definition dtau^2 = -ds^2 is not printed in any local source (geodesic lane doubt). |
| `r_{peri}` | periastron (minimum Schwarzschild radial coordinate) of an unbound timelike geodesic: the largest simple root of R(r) with dR/dr > 0 | cm | `r_p` (gr.timelike-periastron-cubic, gr.timelike-angular-momentum-from-periastron, gr.timelike-proper-time-quadrature; PRD "Periastron distance") | gr.timelike-periastron-cubic | (M1 2026-10-23) Not `r_p`: subscript p already means proton (`m_p`, `n_p`) and polar (`B_p`). |
| `r_{cam}` | Schwarzschild radial coordinate of the camera | cm | `r_{cam}` (gr.bending-integral-finite-radius); `r_c` (pipeline.null-geodesic-trace) | gr.bending-integral-finite-radius | (M1 2026-10-23) Reserved in "Pipeline-level spellings"; now a row. Distinct from `r_0` (Pech83 static observer at infinity) and from `D`. |
| `r_{tp}` | radial turning point (periastron) of a photon orbit with impact parameter b: largest root of r^3 - b^2 r + b^2 r_g = 0 | cm | `r_{tp}` (gr.bending-integral-finite-radius); `u_*` (DA09 Sect. 3, as u = 1/r with M = 1) | gr.bending-integral-finite-radius | (M1 2026-10-23) Exists outside R only for b > b_max (R > 3GM/c^2). Photon analogue of `r_{peri}`. |
| `v^\infty` | asymptotic speed of an unbound timelike geodesic relative to static observers at infinity, c (1 - tilde-E^-2)^(1/2) | cm s^-1 | `v^\infty` (gr.timelike-speed-at-infinity; none printed, derived) | gr.timelike-speed-at-infinity | (M1 2026-10-23) Superscript infinity per the observer-frame convention. |
| `\phi_{Sch}` | Schwarzschild (Boyer-Lindquist at a = 0) azimuthal coordinate | rad | `\phi` (gr.schwarzschild-line-element, gr.kerr-line-element-bl, gr.timelike-constants-schwarzschild, gr.timelike-dphi-dtau-schwarzschild; DA09 eqs. 1, 9) | gr.schwarzschild-line-element | (M1 2026-10-23) Bare `\phi` is the Bel02 surface-element azimuth about the line of sight. For an equatorial orbit `\phi_{Sch}` = `\phi_{orb}` up to the origin. Pairs with `\theta_{Sch}`. |
| `\mathcal{R}(r)` | radial function of the separated geodesic equations (DA09 eq. 6); for a = 0, delta_1 = 1, Q = 0: r^4 E^2 - (r^2 - 2Mr)(L^2 + r^2) = r^4 (dr/dtau)^2 | cm^4 (geometric) | `R(r)` (gr.timelike-radial-equation-schwarzschild, gr.timelike-proper-time-quadrature; DA09 eq. 6) | gr.timelike-radial-equation-schwarzschild | (M1 2026-10-23) Bare `R` is the stellar radius. |
| `\delta_1` | geodesic type switch: 0 null, 1 timelike (DA09) | 1 | `\delta_1` (gr.timelike-constants-schwarzschild, gr.timelike-dt-dtau-schwarzschild, gr.timelike-dphi-dtau-schwarzschild; DA09 eqs. 6-7) | gr.timelike-constants-schwarzschild | (M1 2026-10-23) Distinct from `\delta_{L}` (Reinhard offset) and `\delta_{Dec}` (declination). |
| `a_{spin}` | Kerr parameter; DA09: the angular momentum of the black hole in G = c = M = 1 | cm (geometric) | `a` (gr.kerr-line-element-bl; DA09 eqs. 1-3) | gr.kerr-line-element-bl | (M1 2026-10-23) Kerr parent only, outside PRD scope; set to 0 everywhere in this project. Bare a would collide with `a_{key}`, `a_q`, `a_0`, `a^{B.2}`, `a_{dip}`. |
| `\Delta_K` | Kerr function r^2 - 2Mr + a^2 (DA09 eq. 2 with M = 1) | cm^2 | `\Delta` (gr.kerr-line-element-bl; DA09 eq. 2) | gr.kerr-line-element-bl | (M1 2026-10-23) a = 0: r^2 - r_g r. Kept apart from `\Delta_\rho`, `\Delta T`. |
| `\rho_K^2` | r^2 + a^2 cos^2 theta_Sch (DA09 eq. 2) | cm^2 | `\rho^2` (gr.kerr-line-element-bl; DA09 eq. 2) | gr.kerr-line-element-bl | (M1 2026-10-23) Bare `\rho` is mass density. a = 0: r^2. |
| `\Sigma_K^2` | (r^2 + a^2)^2 - a^2 Delta_K sin^2 theta_Sch (DA09 eq. 3) | cm^4 | `\Sigma^2` (gr.kerr-line-element-bl; DA09 eq. 3) | gr.kerr-line-element-bl | (M1 2026-10-23) a = 0: r^4. |
| `\Theta_K(\theta_{Sch})` | polar function of the separated geodesic equations (DA09 eq. 7), a function of the Boyer-Lindquist polar coordinate theta_Sch | cm^2 (geometric) | `\Theta` (gr.kerr-line-element-bl; DA09 eq. 7) | gr.kerr-line-element-bl | (M1 2026-10-23) Renamed from the proposed `\Theta_K(\theta)`: the argument is the coordinate `\theta_{Sch}`, and bare `\theta` is the RT ray angle. Avoids the old pipeline `Theta` (now `\theta_{Bn}`). |
| `\mathcal{Q}` | Carter constant (DA09) | cm^2 (geometric) | `\mathcal{Q}` (gr.timelike-constants-schwarzschild; DA09 eqs. 6-7) | gr.timelike-constants-schwarzschild | (M1 2026-10-23) 0 for the equatorial orbit chosen at a = 0. Distinct from `Q_\nu` (Stokes) and `Q_{s}` (cohesive energy). |

### Camera frame and image plane (M1 2026-10-23)

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `e^{st}_{(t)}` | time leg of the static orthonormal Schwarzschild tetrad, e = [c (1 - r_g/r)^(1/2)]^(-1) d/dt (= u^mu_obs/c) | s cm^-1 (component along d/dt) | `e_(t)` (gr.static-tetrad; MG09 eq. 2.2.6); FIDO time leg (James+15, a = 0) | gr.static-tetrad | (M1 2026-10-23) Label st = local static observer, the same observer as `u^\mu_{obs}` (Pech83). |
| `e^{st}_{(r)}` | radial leg of the static tetrad, (1 - r_g/r)^(1/2) d/dr | 1 | `e_(r)` (gr.static-tetrad); `e_r-hat` (James+15 eq. A.3, gr.camera-aberration); `e_lambda-hat` (WS83, gr.static-tetrad-field-projection) | gr.static-tetrad | (M1 2026-10-23) |
| `e^{st}_{(\theta)}` | polar leg of the static tetrad, r^(-1) d/dtheta_Sch | cm^-1 | `e_(vartheta)` (gr.static-tetrad); `e_theta-hat` (James+15 eq. A.3) | gr.static-tetrad | (M1 2026-10-23) Leg label refers to `\theta_{Sch}`, not the RT ray angle `\theta`. |
| `e^{st}_{(\phi)}` | azimuthal leg of the static tetrad, (r sin theta_Sch)^(-1) d/dphi_Sch | cm^-1 | `e_(varphi)` (gr.static-tetrad); `e_phi-hat` (James+15 eq. A.3) | gr.static-tetrad | (M1 2026-10-23) Leg label refers to `\phi_{Sch}`, not the surface-element azimuth `\phi`. |
| `e_{(a)}` | camera orthonormal tetrad leg, a in {0, 1, 2, 3}; e_(0) = u_cam/c | e_(0)^t in s cm^-1; spatial legs as the static legs they reduce to | `e_(i)`, `e-hat_(i)` (gr.tetrad-lorentz-transformation; MG09 eq. 1.4.6); `e_x, e_y, e_z` (James+15 Sect. 2.1, gr.camera-aberration); `e_(a)` (pipeline.camera-tetrad-projection) | gr.tetrad-lorentz-transformation | (M1 2026-10-23) Reserved in "Pipeline-level spellings"; now a row. Not the polarisation vector `\mathbf{e}_j` or the charge `e`. |
| `e^{(a)}` | dual (co)tetrad, e^(a)_mu e_(b)^mu = delta^a_b | inverse of e_(a) | `theta^(i)` (gr.tetrad-vector-components; MG09 eqs. 1.4.1, 1.4.19) | gr.tetrad-vector-components | (M1 2026-10-23) MG09 theta^(i) renamed: `\theta` is the RT ray angle. |
| `u^\mu_{cam}` | four-velocity of the camera, u = c gamma_cam (e_(0)^st + beta_cam v-hat_cam) | as `u^\mu_{obs}` (u^t dimensionless with t in s; u^r cm s^-1; u^theta, u^phi s^-1) | `u` (gr.camera-four-velocity-static-frame; MG09 eq. 1.4.18); `u^mu` (pipeline.camera-worldline) | gr.camera-four-velocity-static-frame | (M1 2026-10-23) Reserved in "Pipeline-level spellings"; now a row. |
| `\beta_{cam}` | camera speed relative to the local static observer, in units of c | 1 | `beta` (gr.camera-four-velocity-static-frame, gr.tetrad-lorentz-transformation; MG09 eq. 1.4.7; James+15 eq. A.9, speed relative to the FIDO) | gr.camera-four-velocity-static-frame | (M1 2026-10-23) Not the bending angle `\beta`. |
| `\gamma_{cam}` | Lorentz factor of the camera relative to the local static observer, (1 - beta_cam^2)^(-1/2) | 1 | `gamma` (gr.tetrad-lorentz-transformation; MG09 eq. 1.4.7) | gr.tetrad-lorentz-transformation | (M1 2026-10-23) Not `\gamma_B` or `\gamma^2`. |
| `\hat v_{cam}` | unit 3-vector of the camera's direction of motion in the local static frame | 1 | `n^a`, `n` (gr.tetrad-lorentz-transformation, gr.camera-four-velocity-static-frame; MG09 eqs. 1.4.7, 1.4.18); `B` (James+15 Sect. 2.1, e_y = B) | gr.tetrad-lorentz-transformation | (M1 2026-10-23) MG09 n renamed (bare n retired); James's B is not a magnetic field. |
| `\Lambda_{boost}` | Lorentz boost matrix between the static and camera tetrads | 1 | `Lambda^a_b`, `A_i^k` (gr.tetrad-lorentz-transformation; MG09 eqs. 1.4.6-1.4.7) | gr.tetrad-lorentz-transformation | (M1 2026-10-23) Not the metric function `\Lambda`. Sign convention: see gr.tetrad-lorentz-transformation notes. |
| `\hat N_{cam}` | unit vector in the camera frame pointing from the camera toward the apparent source of an incoming ray (opposite to propagation) | 1 | `N`, `N_x, N_y, N_z` (gr.camera-aberration; James+15 eq. A.8) | gr.camera-aberration | (M1 2026-10-23) Photon propagation direction in the camera frame is -N_cam. Not `N_A`, `N_H`, `N_{pix}`. |
| `\hat n_{st}` | propagation direction of an incoming ray measured by the local static observer | 1 | `n_F` (gr.camera-aberration; James+15 eq. A.9) | gr.camera-aberration | (M1 2026-10-23) Decorated n (hat and subscript); bare n stays retired. |
| `\theta_{cam}` | polar angle on the camera's local sky, measured from camera axis e_(3) (James e_z) | rad | `theta_cs` (gr.camera-aberration; James+15 eq. A.8) | gr.camera-aberration | (M1 2026-10-23) |
| `\phi_{cam}` | azimuth on the camera's local sky about e_(3), from e_(1) (James e_x) | rad | `phi_cs` (gr.camera-aberration; James+15 eq. A.8) | gr.camera-aberration | (M1 2026-10-23) |
| `\nu_{cam}` | photon frequency measured in the camera frame | Hz | `f_c` (gr.camera-blueshift-celestial, gr.doppler-factor-sr; James+15 eq. A.16); `nu_0` (gr.frequency-ratio-general; DA09 eq. 70, observed) | gr.frequency-ratio-general | (M1 2026-10-23) Reserved in "Pipeline-level spellings"; now a row. DA09 nu_0 is the observed frequency, Bel02 nu_0 the emitter frequency (`\nu`): the source spelling nu_0 has two meanings. |
| `\nu_{st}` | photon frequency measured by the local static observer at the camera's position | Hz | text/other: FIDO frequency (gr.doppler-factor-sr; James+15 eq. A.11, E_F) | gr.doppler-factor-sr | (M1 2026-10-23) |
| `g_{shift}` | frequency ratio nu_cam/nu = (p.u_cam)/(p.u_em), with nu the emitter-frame frequency | 1 | `g` (gr.frequency-ratio-general, gr.intensity-g-cubed; DA09 eq. 70, g = nu_0/nu); `f_c/f'` (James+15 eq. A.16, emitter at rest at infinity) | gr.frequency-ratio-general | (M1 2026-10-23) Reserved in "Pipeline-level spellings"; bare g is surface gravity. Static emitter at R, camera at infinity: g_shift = (1 - x_g)^(1/2) = 1/(1 + z_g). |
| `\mathcal{D}_{cam}` | special-relativistic Doppler factor nu_cam/nu_st = (1 - beta_cam^2)^(1/2)/(1 - beta_cam N_y) | 1 | first factor of `f_c/f'` (gr.doppler-factor-sr, gr.camera-blueshift-celestial; James+15 eq. A.16) | gr.doppler-factor-sr | (M1 2026-10-23) `D` is distance. |
| `f_{cam}` | camera focal length in pixel units (projection constant) | 1 (pixels) | none (project definition; gr.camera-projection-convention) | gr.camera-projection-convention | (M1 2026-10-23) Requires f_cam > 0. Not `f_\nu`, `f_{col}`, `f_B`. The physical focal length is `f_{foc}` (cm); f_cam = f_foc / pixel pitch (pixel pitch has no row yet). |
| `r_{pix}` | radial distance of a pixel from the optical-axis pixel | 1 (pixels) | none (project definition; gr.camera-projection-convention) | gr.camera-projection-convention | (M1 2026-10-23) Pixel label as in `x_{pix}`, `y_{pix}`. |
| `\theta_{field}` | field angle between N_cam and the optical axis a_cam | rad | none (project definition; gr.camera-projection-convention) | gr.camera-projection-convention | (M1 2026-10-23) Pinhole: tan theta_field = r_pix/f_cam; equidistant: theta_field = r_pix/f_cam (open PRD decision). Requires f_cam > 0 and 0 <= theta_field <= pi; the full field of view may exceed pi. Not the PSF offset `\theta_{psf}`. |
| `\phi_{field}` | image-plane azimuth of a pixel about the optical axis, from image-right | rad | none (project definition; gr.camera-projection-convention) | gr.camera-projection-convention | (M1 2026-10-23) |
| `\hat a_{cam}` | optical-axis (pointing) unit vector in the camera frame | 1 | none (project definition; PRD The shot, Pointing; gr.camera-projection-convention) | gr.camera-projection-convention | (M1 2026-10-23) Visually close to the vacuum coefficient `\hat a`; distinct string. |

### EOS fits (BSk)

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `\xi_\rho` | xi = log10(rho / g cm^-3); log denotes log10, ln the natural logarithm | 1 | `xi` (eos.notation-bsk-fit) | eos.notation-bsk-fit | log = log10, ln = natural log (BSk convention, not applied elsewhere). |
| `\zeta_P` | zeta = log10(P / dyn cm^-2) | 1 | `zeta` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `n_b` | baryon number density | cm^-3 | `n_b` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `n_{fm}` | dimensionless baryon density n = n_b / fm^-3 used in the fits | 1 | `n` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `\rho_s` | value of rho at the neutron-star surface, set to the density of 56Fe at zero pressure and temperature | g cm^-3 | `rho_s` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `n_s` | value of n_b at the neutron-star surface | cm^-3 | `n_s` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `m_u^{BSk}` | atomic mass unit as used in the fits, 1.66e-24 g | g | `m_u` (eos.notation-bsk-fit) | eos.notation-bsk-fit | Fit-internal rounded value 1.66e-24 g only. |
| `\mathcal{H}` | pseudo-enthalpy H = int_0^P dP'/(rho(P') c^2 + P') = ln(h/h_s) | 1 | `H` (eos.notation-bsk-fit) | eos.notation-bsk-fit | BSk pseudo-enthalpy; bare H is not used (H_nu is the Eddington flux). |
| `h_b` | zero-temperature enthalpy per baryon h = (rho c^2 + P)/n_b | erg | `h` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `h_{b,s}` | value of h at the stellar surface | erg | `h_s` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `\eta_{\mathcal{H}}` | eta = exp(H) - 1 | 1 | `eta` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `\Delta_\rho` | fractional deviation from the linear law rho = n_b m_u, Delta_rho = rho/(n_b m_u) - 1 | 1 | `Delta_rho` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `\tilde\rho` | dimensionless density rho~ = (rho/m_u) fm^3 = rho/(1.66e15 g cm^-3) | 1 | `rho~` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `a_i^{BSk}` | fit parameters of Eq. (3), i = 1..23, per EoS model (Table 2) | 1 | `a_i` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `b_i^{BSk}` | fit parameters of Eq. (6), i = 1..8, per EoS model (Table 3) | 1 | `b_i` (eos.notation-bsk-fit) | eos.notation-bsk-fit |  |
| `c_i^{BSk}` | fit parameters of Eq. (7), i = 1..9, per EoS model (Table 4) | 1 | `c_i` (eos.notation-bsk-fit) | eos.notation-bsk-fit | Not the radiation constants c_1, c_2. |
| `\Gamma_{ad}` | adiabatic index Gamma = (n_b/P) dP/dn_b (note: unrelated to the plasma Coulomb coupling parameter) | 1 | `Gamma` (eos.notation-bsk-fit) | eos.notation-bsk-fit | Not the Coulomb coupling Gamma_coul. |

### Plasma, magnetic field and cyclotron quantities

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `E` | photon energy hbar omega | keV (code: erg) | `E_{ph}` (rt.notation-cond-energies-geometry); `E` (rt.notation-val) | rt.notation-val |  |
| `B_{14}` | local magnetic field strength B/(10^14 G) | 1 | `B_{14}` (rt.notation-cond-energies-geometry); `B_14` (opacity.notation-fields, rt.notation-vacuum) | opacity.notation-fields | (r1 fix applied: B_14 = B/(1e14 G) taken from van_adelsberg_lai_2006_magnetic_atm.tex L360 "B_{14}=B/(10^{14}\ \mbox{G})" (txt L185 "B14 = B/(10 G), and", exponent lost by the text layer); recorded in opacity.notation-fields notes, truth/entries/opacity.toml:33; entry source stays PC03, no per-symbol source field.) |
| `E_{ce}` | electron cyclotron energy hbar omega_ce = hbar e B/(m_e c) = 1158 B_14 keV | keV | `E_{ce}` (rt.notation-cond-energies-geometry); `E_Be` (rt.notation-vacuum); text/other: `E_{Be}` | rt.notation-vacuum |  |
| `E_{ci}` | ion cyclotron energy hbar omega_ci = hbar Z e B/(A m_u c) = 0.635 (Z/A) B_14 keV | keV | `E_{ci}` (rt.notation-cond-energies-geometry); `E_ci` (rt.notation-spw09); `E_Bi` (rt.notation-vacuum); text/other: `E_{Bi}` | rt.notation-spw09 |  |
| `A_{ion}` | ion mass number A | 1 | `A` (rt.notation-cond-matter-dielectric); `A_ion` (opacity.notation-particles-plasma); text/other: `A_{mass}` | rt.notation-cond-matter-dielectric | (r1 fix applied: vAL05 spells it A; A_mass is a historical project alias. Approximate atomic-mass usage, including noninteger examples, as in the source.) |
| `Z_{ion}` | ion charge number Z (nuclear charge of a fully ionised ion) | 1 | `Z` (rt.notation-cond-matter-dielectric); `Z_ion` (opacity.notation-particles-plasma); text/other: `Z_{nuc}` | rt.notation-cond-matter-dielectric | Bare Z is the CIE tristimulus value. (r1 fix applied: vAL05 spells it Z; Z_nuc is a historical project alias.) |
| `n_e` | electron number density | cm^-3 | `n_{e}` (rt.notation-cond-matter-dielectric); `n_e` (opacity.notation-particles-plasma) | rt.notation-cond-matter-dielectric |  |
| `b_{rel}` | relativistic field parameter b = B/B_Q = hbar omega_ce/(m_e c^2) = B_12/44.14 | 1 | `b` (eos.notation-electron-gas); `b_rel` (opacity.notation-fields) | eos.notation-electron-gas |  |
| `\omega_{ce}` | electron cyclotron angular frequency omega_ce = eB/(m_e c) (Gaussian CGS) | s^-1 | `omega_c` (eos.notation-electron-gas); `omega_ce` (opacity.notation-photon-cyclotron) | eos.notation-electron-gas |  |
| `B_{12}` | magnetic field B/(10^12 G) | 1 | `B_12` (eos.notation-electron-gas, opacity.notation-fields) | eos.notation-electron-gas |  |
| `\lambda_e` | electron Compton wavelength hbar/(m_e c) | cm | `lambda_e` (eos.notation-electron-gas) | eos.notation-electron-gas |  |
| `a_m` | magnetic length (hbar c/(eB))^(1/2) = lambda_e / sqrt(b) | cm | `a_m` (eos.notation-electron-gas, opacity.notation-fields) | eos.notation-electron-gas |  |
| `\mu_e` | electron chemical potential including rest energy | erg | `mu_e (source mu)` (eos.notation-electron-gas); text/other: `mu` | eos.notation-electron-gas | Electron chemical potential. Not the ray cosine mu, not mean molecular weight. |
| `n_L` | Landau level index n = 0..n_max in the electron-pressure sum | 1 | `n (Landau)` (eos.notation-electron-gas) | eos.notation-electron-gas |  |
| `\chi_n` | chi_n = (mu - m_e c^2 sqrt(1+2bn))/(k_B T), reduced chemical potential of level n | 1 | `chi_n` (eos.notation-electron-gas) | eos.notation-electron-gas |  |
| `\tau_n` | tau_n = (1/sqrt(1+2bn)) T/T_r, reduced temperature of level n | 1 | `tau_n` (eos.notation-electron-gas) | eos.notation-electron-gas |  |
| `P_r` | relativistic unit of pressure m_e c^2/lambda_e^3 = 1.4218e25 dyn cm^-2 | dyn cm^-2 | `P_r` (eos.notation-electron-gas) | eos.notation-electron-gas |  |
| `T_r` | relativistic unit of temperature m_e c^2/k_B = 5.930e9 K | K | `T_r` (eos.notation-electron-gas) | eos.notation-electron-gas |  |
| `I_{1/2}(\chi,\tau)` | standard Fermi-Dirac integral of order 1/2 for relativistic electrons (given by a fit in Chabrier & Potekhin 1998 per PY01) | 1 | `I_{1/2}(chi, tau)` (eos.notation-electron-gas) | eos.notation-electron-gas |  |
| `B` | magnetic field strength | G | `B` (opacity.notation-fields) | opacity.notation-fields |  |
| `B_0` | atomic magnetic field m_e^2 c e^3 / hbar^3 = 2.3505e9 G | G | `B_0` (opacity.notation-fields) | opacity.notation-fields |  |
| `B_Q` | QED critical field m_e^2 c^3 / (e hbar) = 4.414e13 G (called B_r in PC03) | G | `B_Q` (opacity.notation-fields) | opacity.notation-fields |  |
| `\gamma_B` | atomic field parameter B / B_0 (PHC16 writes gamma) | 1 | `gamma_B` (opacity.notation-fields) | opacity.notation-fields | (r1 fix applied: gamma = B/B_0 from potekhin_ho_chabrier_2016_review.txt L120-121 "γ ≡ B/B0 ." (source spelling gamma); recorded in opacity.notation-fields notes, truth/entries/opacity.toml:33; no per-symbol source field.) |
| `\mu_{red}` | electron-proton reduced mass m_e m_p / m_H (PC03 writes mu) | g | `mu_red` (opacity.notation-particles-plasma); text/other: `mu` | opacity.notation-particles-plasma |  |
| `Y_e` | electron fraction Z/A | 1 | `Y_e` (opacity.notation-particles-plasma); `Y_e (VAL06)` (rt.notation-vacuum) | rt.notation-vacuum |  |
| `T` | gas temperature | K | `T` (opacity.notation-particles-plasma) | opacity.notation-particles-plasma |  |
| `T_6` | T / (10^6 K) | 1 | `T_6` (opacity.notation-particles-plasma, rt.notation-vacuum) | opacity.notation-particles-plasma |  |
| `\rho_1` | rho/(1 g cm^-3), dimensionless density | 1 | `rho_0` (opacity.notation-particles-plasma); `ρ_1` (rt.notation-vacuum) | opacity.notation-particles-plasma |  |
| `n_p` | proton number density | cm^-3 | `n_p` (opacity.notation-particles-plasma) | opacity.notation-particles-plasma |  |
| `x_H` | number fraction of optically active bound hydrogen atoms in PC03 eq. 31 (0 for a fully ionised plasma); dissolved and highly excited states are excluded, so it is not silently the thermodynamic neutral fraction | 1 | `x_H` (opacity.notation-particles-plasma) | opacity.notation-particles-plasma |  |
| `f_ion` | ionisation fraction (McPHAC notation) | 1 | `f_ion` (opacity.notation-particles-plasma) | opacity.notation-particles-plasma | (r1 fix applied: McPHAC definition at haakonsen_2012_mcphac.txt L225 "For an ionization fraction fion at a given depth, the reduced" (eq. 11, L225-232); recorded in opacity.notation-particles-plasma notes, truth/entries/opacity.toml:74; no per-symbol source field.) |
| `\omega_{pe}` | electron plasma angular frequency (4 pi e^2 n_e / m_e)^(1/2) (PHC16 writes omega_pe) | rad/s | `omega_pl` (opacity.notation-particles-plasma) | opacity.notation-particles-plasma |  |
| `\rho_B` | density below which a strongly quantising field keeps all electrons in the ground Landau level (requires also zeta_e >> 1, i.e. low enough temperature) | g/cm^3 | `rho_B` (opacity.notation-particles-plasma) | opacity.notation-particles-plasma |  |
| `\omega` | photon angular frequency | rad/s | `omega` (opacity.notation-photon-cyclotron) | opacity.notation-photon-cyclotron |  |
| `u_{kT}` | hbar omega / (k_B T), scaled photon energy (paper symbol u; McPHAC: u = h nu / k_B T) | 1 | `u_kT` (opacity.notation-photon-cyclotron); text/other: `u` | opacity.notation-photon-cyclotron | Bare u is retired (atomic mass unit, compactness, Pechenick M/r, hbar omega/kT all collided). (r1 fix applied: opacity.notation-photon-cyclotron locator now 'eq. 1; eq. 41 (u = hbar omega/(k_B T)); ...'; potekhin_chabrier_2003_ff_opacity.txt L659-664 "1 − e−u Λffα , u ≡"; truth/entries/opacity.toml:102.) |
| `\omega_{cp}` | proton cyclotron angular frequency (m_e/m_p) omega_ce | rad/s | `omega_cp` (opacity.notation-photon-cyclotron) | opacity.notation-photon-cyclotron |  |
| `\omega_{ci}` | ion cyclotron angular frequency Z e B / (m_i c) | rad/s | `omega_ci` (opacity.notation-photon-cyclotron) | opacity.notation-photon-cyclotron | (r1 fix applied: potekhin_ho_chabrier_2016_review.txt L332-333 "ZeB/(mi c) is the ion cyclotron frequency"; recorded in notes of truth/entries/opacity.toml:102; no per-symbol source field.) |
| `\zeta_e` | electron quantisation parameter hbar omega_ce/(k_B T) | 1 | `beta_e` (opacity.notation-photon-cyclotron); text/other: `beta_e (PC03)`, `zeta_e (PHC16)` | opacity.notation-photon-cyclotron |  |
| `\zeta_p` | hbar omega_cp / (k_B T), proton quantisation parameter | 1 | `beta_p` (opacity.notation-photon-cyclotron) | opacity.notation-photon-cyclotron | (r1 fix applied: locator now cites the proton-quantisation discussion, potekhin_chabrier_2003_ff_opacity.txt L715-716 "βp < 1, and obtain"; truth/entries/opacity.toml:102.) |
| `\zeta_i` | hbar omega_ci / (k_B T), ion quantisation parameter (PHC16) | 1 | `zeta_i` (opacity.notation-photon-cyclotron) | opacity.notation-photon-cyclotron | (r1 fix applied: potekhin_ho_chabrier_2016_review.txt L332-333 "ζi = h̄ωci /kB T = 0.0737 (Z/A)B12 /T6"; recorded in notes of truth/entries/opacity.toml:102; no per-symbol source field.) |

### Normal modes, dielectric tensor and vacuum polarisation

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `\varepsilon_\perp` | transverse diagonal element of the dielectric tensor I + 4 pi chi, frame z || B: eps_xx = eps_yy = eps_perp | 1 | `ε` (rt.notation-cond-matter-dielectric, opacity.notation-modes-polarization); text/other: `epsilon` | rt.notation-cond-matter-dielectric |  |
| `\varepsilon_\wedge` | off-diagonal (gyrotropic) element of the dielectric tensor I + 4 pi chi, frame z || B, with eps_xy = +i eps_wedge and eps_yx = -i eps_wedge (convention common to vAL05 and PLCH04) | 1 | `g` (rt.notation-cond-matter-dielectric, opacity.notation-modes-polarization) | rt.notation-cond-matter-dielectric |  |
| `\varepsilon_\parallel` | longitudinal element (along B) of the dielectric tensor I + 4 pi chi, frame z || B: eps_zz = eps_parallel | 1 | `η` (rt.notation-cond-matter-dielectric, opacity.notation-modes-polarization); text/other: `eta` | rt.notation-cond-matter-dielectric |  |
| `\theta_{kB}` | angle between the photon wave vector k and the magnetic field B | rad | `theta_B` (opacity.notation-modes-polarization); `α_Bk` (rt.notation-spw09); `θ_kB` (rt.notation-val) | rt.notation-val | PC03/PLCH04 write theta_B for this; renamed so theta_B is the magnetic colatitude. |
| `\alpha_{pol}` | basic polarisation index, 0 (linear along B) or +/-1 (circular in the transverse plane); PC03 and PLCH04 write alpha | 1 | `alpha_pol` (opacity.notation-modes-polarization); text/other: `alpha` | opacity.notation-modes-polarization | Basic (cyclic) polarisation index -1, 0, +1. Inside subscripts of opacity quantities (sigma_alpha, kappa-hat_alpha, nu_alpha) the dummy index alpha means alpha_pol. |
| `j` | normal-mode index, j = 1 extraordinary (X), j = 2 ordinary (O) | 1 | `j` (opacity.notation-modes-polarization); `i` (rt.notation-spw09) | opacity.notation-modes-polarization | Normal-mode index j = 1 extraordinary (X), j = 2 ordinary (O); i, j for pairs. Never an angle-node index (q) or a linear-polarisation index (l). (r1 fix applied: owner column changed to opacity.notation-modes-polarization; evidence potekhin_lai_chabrier_ho_2004_polarization.txt L355 "where j = 1, 2 correspond to the extraordinary mode (X-", recorded in truth/entries/opacity.toml:138.) |
| `\mathbf{e}_j` | polarisation vector of normal mode j (unit complex vector) | 1 | `e_j` (opacity.notation-modes-polarization) | opacity.notation-modes-polarization |  |
| `e_{j,\alpha}` | cyclic component of e_j: e_{j,0} = e_{j,z}, e_{j,+/-1} = (e_{j,x} +/- i e_{j,y})/sqrt(2) | 1 | `e_j_alpha` (opacity.notation-modes-polarization) | opacity.notation-modes-polarization | (r1 fix applied: PC03 potekhin_chabrier_2003_ff_opacity.txt L382-383 "ej,±1 = (ej,x ±iej,y )/ 2 are the circular components." added to notes of truth/entries/opacity.toml:138; no per-symbol source field.) |
| `K_j` | transverse polarisation parameter of mode j in the frame with z along k and B in the x-z plane | 1 | `K_j` (opacity.notation-modes-polarization) | opacity.notation-modes-polarization |  |
| `K_{z,j}` | longitudinal polarisation parameter K_{z,j} of mode j | 1 | `K_z_j` (opacity.notation-modes-polarization) | opacity.notation-modes-polarization |  |
| `\beta_{pol}` | mode polarisation parameter beta of PLCH04 Eq. 21 (not the quantisation parameter beta_e) | 1 | `beta_pol` (opacity.notation-modes-polarization); text/other: `beta` | opacity.notation-modes-polarization |  |
| `\varepsilon'_\perp` | epsilon' = eps_perp + a_hat, transverse element including vacuum | 1 | `eps_prime` (opacity.notation-modes-polarization) | opacity.notation-modes-polarization |  |
| `\varepsilon'_\parallel` | eta' = eta_par + a_hat + q_vac, longitudinal element including vacuum | 1 | `eta_prime` (opacity.notation-modes-polarization) | opacity.notation-modes-polarization |  |
| `\chi` | complex polarisability tensor of the plasma, chi = chi^H + i chi^A | 1 | `chi` (opacity.notation-modes-polarization) | opacity.notation-modes-polarization |  |
| `\chi^H` | Hermitian part of chi (refraction and polarisation) | 1 | `chi_H` (opacity.notation-modes-polarization) | opacity.notation-modes-polarization |  |
| `\chi^A` | chi^A in chi = chi^H + i chi^A; i chi^A is the anti-Hermitian part (absorption) | 1 | `chi_A` (opacity.notation-modes-polarization) | opacity.notation-modes-polarization |  |
| `\chi^{vac}` | polarisability tensor of the vacuum | 1 | `chi_vac` (opacity.notation-modes-polarization) | opacity.notation-modes-polarization |  |
| `\hat a` | vacuum polarisation coefficient a-hat (isotropic part of 4 pi chi_vac); vanishes at B = 0 | 1 | `a_hat` (opacity.notation-vacuum) | opacity.notation-vacuum |  |
| `q_{vac}` | vacuum polarisation coefficient q of 4 pi chi_vac = diag(a_hat, a_hat, a_hat + q) (PLCH04 Eq. 3) | 1 | `q_vac` (opacity.notation-vacuum) | opacity.notation-vacuum | (r1 fix applied: locator of opacity.notation-vacuum now includes PLCH04 eq. 3, potekhin_lai_chabrier_ho_2004_polarization.txt L163 "4πχvac = diag(â, â, â + q)"; truth/entries/opacity.toml:162.) |
| `m_{vac}` | vacuum magnetisation coefficient m of mu^-1 = I + diag(a_hat, a_hat, a_hat + m) (PLCH04 Eq. 4) | 1 | `m_vac` (opacity.notation-vacuum) | opacity.notation-vacuum | (r1 fix applied: locator of opacity.notation-vacuum now includes PLCH04 eq. 4, potekhin_lai_chabrier_ho_2004_polarization.txt L78 "µ−1 = I + diag(â, â, â + m)"; truth/entries/opacity.toml:162.) |
| `\xi_{vac}` | inverse relativistic field parameter 1/b (paper symbol xi) | 1 | `xi_vac` (opacity.notation-vacuum); text/other: `xi` | opacity.notation-vacuum |  |
| `E_V` | photon energy of the vacuum resonance (plasma and vacuum contributions to the dielectric tensor cancel) at given B and rho | keV | `E_V(ρ)` (rt.notation-spw09); `E_V (VAL06)` (rt.notation-vacuum) | rt.notation-vacuum | (r1 fix: the PLCH04 declaration (opacity.notation-vacuum) was removed; the definition is VAL06 (rt.notation-vacuum).) |
| `\rho_V` | density at which a photon of energy E meets the vacuum resonance | g cm^-3 | `ρ_V` (rt.notation-vacuum); text/other: `rho_V` | rt.notation-vacuum | (r1 fix: the PLCH04 declaration (opacity.notation-vacuum) was removed; the definition is VAL06 (rt.notation-vacuum).) |
| `f_B` | slowly varying function of B, f(B) ~ [3 delta_V/(q + m)]^(1/2) ~ 1 (Ho & Lai 2003, lines 444-457 of the local copy; VAL06 uses f_B ~ 1) | 1 | `f_B` (opacity.notation-vacuum); `f_B (VAL06)` (rt.notation-vacuum) | rt.notation-vacuum | Defined only in Ho & Lai 2003 eq. 2.41 (not a licence to set f_B = 1). |
| `E_1` | photon energy E / (1 keV) | 1 | `E_1` (rt.notation-vacuum) | opacity.notation-vacuum | (r1 fix: the PLCH04 declaration was removed; E_1 = E/(1 keV) is defined in VAL06 (rt.notation-vacuum).) |
| `E_{ad}` | adiabatic energy; resonant conversion occurs when E >~ E_ad | keV | `E_ad` (rt.notation-vacuum) | rt.notation-vacuum |  |
| `P_{jump}` | non-adiabatic jump probability; 1-P_jump is the O<->X conversion probability | 1 | `P_jump` (rt.notation-vacuum) | rt.notation-vacuum |  |
| `H_\rho` | density scale height along the ray, \|ds/d ln ρ\| | cm | `H_ρ` (rt.notation-vacuum) | rt.notation-vacuum |  |
| `B_l` | approximate, energy- and temperature-dependent field scale above which the vacuum resonance lies between the X- and O-mode photospheres for typical photons (not a universal field threshold) | G | `B_l` (rt.notation-vacuum) | rt.notation-vacuum |  |
| `S_{stim}` | stimulated-emission factor S_stim = 1 - exp(-E/k_B T) | 1 | `S` (rt.notation-vacuum) | rt.notation-vacuum |  |

### Opacities and cross sections

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `\kappa_j` | total extinction opacity of normal mode j, kappa^a_j + kappa^s_j | cm^2 g^-1 | `kappa_j` (opacity.notation-opacity-symbols); `κ_j^tot` (rt.notation-val) | rt.notation-val |  |
| `\kappa^a_j` | true-absorption (free-free for fully ionised H) opacity of normal mode j; depends on direction through theta_kB | cm^2 g^-1 | `kappa_a_j` (opacity.notation-opacity-symbols); `k_ν^i` (rt.notation-spw09); `κ_j^abs` (rt.notation-val); text/other: `kappa^abs_j`, `κ_j^ff` | rt.notation-val |  |
| `\kappa^s_{jj'}` | scattering opacity kappa^s_{jj'} from mode j into mode j' | cm^2/g | `kappa_s_jj` (opacity.notation-opacity-symbols); `κ_ji^sc` (rt.notation-val); text/other: `kappa^sc_jj'`, `κ_ji^sc (j -> i)` | opacity.notation-opacity-symbols |  |
| `\kappa^s_j` | total scattering opacity out of normal mode j, kappa^s_{j1} + kappa^s_{j2}; depends on direction through theta_kB | cm^2 g^-1 | `kappa_s_j` (opacity.notation-opacity-symbols); `σ_ν^i` (rt.notation-spw09); `κ_j^sc` (rt.notation-val) | rt.notation-spw09 |  |
| `\hat\kappa_\alpha` | basic-polarisation opacity kappa-hat_alpha, mu_alpha = rho kappa-hat_alpha, independent of theta_B (PLCH04 eq. 7) | cm^2/g | `kappa_hat_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols | (r1 fix applied: mu_alpha = rho kappa-hat_alpha, potekhin_lai_chabrier_ho_2004_polarization.txt L114-127 (eq. 7, "α ≡ ρκ̂α (ω) ="); recorded in notes of truth/entries/opacity.toml:206; no per-symbol source field.) |
| `\kappa_\nu` | B = 0 true-absorption opacity at photon frequency nu (free-free for fully ionised hydrogen; McPHAC also allows other absorption processes) | cm^2 g^-1 | `kappa_nu` (opacity.notation-opacity-symbols); `κ_ν` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `\sigma_{T,r}` | reduced Thomson scattering opacity of McPHAC (per gram, includes the ionisation fraction) | cm^2 g^-1 | `sigma_T_r` (opacity.notation-opacity-symbols); `σ_T,r` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `\sigma^a_\alpha` | absorption cross section per hydrogen nucleus for basic polarisation alpha, including the weighted bound and free contributions | cm^2 | `sigma_a_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\sigma^s_\alpha` | scattering cross section for basic polarisation alpha, sigma^{s,e}_alpha + sigma^{s,p}_alpha | cm^2 | `sigma_s_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\sigma^{s,e}_\alpha` | electron scattering cross section sigma^{s,e}_alpha | cm^2 | `sigma_se_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\sigma^{s,p}_\alpha` | proton scattering cross section sigma^{s,p}_alpha | cm^2 | `sigma_sp_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\sigma^{ff}_\alpha` | free-free absorption cross section sigma^{ff}_alpha from electron-proton collisions | cm^2 | `sigma_ff_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\nu^{ff}_\alpha` | effective electron-proton collision frequency for absorption, nu^{ff}_alpha | 1/s | `nu_ff_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\nu_{e,\alpha}` | electron damping frequency nu_{e,alpha} = nu^{ff}_alpha + nu^s_e | 1/s | `nu_e_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\nu^s_e` | natural radiative width of the electron cyclotron resonance, 2 e^2 omega^2/(3 m_e c^3) | 1/s | `nu_s_e` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\nu^s_p` | natural radiative width of the proton cyclotron resonance, 2 e^2 omega^2/(3 m_p c^3) | 1/s | `nu_s_p` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols | (r1 fix applied: locator of opacity.notation-opacity-symbols now includes eqs. 52-53; PC03 L817-822 "νp,α = ναff,p + νps + ναpp"; truth/entries/opacity.toml:206.) |
| `\nu_{p,\alpha}` | proton damping frequency nu_{p,alpha} | 1/s | `nu_p_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols | (r1 fix applied: locator now includes PC03 eq. 53, L817-822 "νp,α = ναff,p + νps + ναpp"; truth/entries/opacity.toml:206.) |
| `\tilde\nu_\alpha` | combined damping frequency nu-tilde_alpha of PC03 eq. 52 | 1/s | `nu_tilde_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols | (r1 fix applied: locator now includes PC03 eq. 52, L750-756 "ν̃α = ναff + (1 + αωce /ω) νp + (1 − αωcp /ω) νe"; truth/entries/opacity.toml:206.) |
| `\Lambda^{ff}_\alpha` | Coulomb logarithm Lambda^{ff}_alpha for free-free absorption (Born, infinite proton mass, fixed Coulomb centre) | 1 | `Lambda_ff_alpha` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\Lambda^{ff}_{cl}` | classical (B = 0, Born) Coulomb logarithm e^{u/2} K_0(u/2) | 1 | `Lambda_ff_cl` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\bar g_{ff}` | temperature-averaged free-free Gaunt factor g-bar_ff; Lambda = (pi/sqrt(3)) g-bar | 1 | `gaunt_ff` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\bar g_{Born}` | Born-approximation Gaunt factor g-bar_Born = (sqrt(3)/pi) Lambda^{ff}_cl | 1 | `gaunt_Born` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols |  |
| `\gamma^2` | scaled temperature gamma^2 = Z^2 Ry / (k_B T) used to index the Gaunt-factor table (McPHAC) | 1 | `gamma_sq` (opacity.notation-opacity-symbols) | opacity.notation-opacity-symbols | (r1 fix applied: McPHAC haakonsen_2012_mcphac.txt L246-249 "and scaled temperature γ 2 = Z 2 Ry/(kB T )"; recorded in notes of truth/entries/opacity.toml:206; no per-symbol source field.) |
| `k_\nu` | total opacity k_ν = κ_ν + σ_T,r (called total absorption opacity in the paper) | cm^2 g^-1 | `k_ν` (rt.notation-mcphac) | rt.notation-mcphac | McPHAC B = 0 total opacity kappa_nu + sigma_T,r. SPW09 k_nu^i (absorption) is renamed kappa^a_j. |
| `\kappa_T` | Thomson scattering opacity, 0.4 cm^2 g^-1 | cm^2 g^-1 | `κ_T` (rt.notation-val) | rt.notation-val | Fixed 0.4 cm^2 g^-1 (VAL06); not sigma_T,r. |

### Condensed surface (vAL05, Pot12)

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `B_{13}` | local magnetic field strength in units of 1e13 G (B/1e13 G) | 1 | `B_{13}` (rt.notation-cond-energies-geometry) | rt.notation-cond-energies-geometry |  |
| `E_{pe}` | electron plasma energy of the condensed matter | keV | `E_{pe}` (rt.notation-cond-energies-geometry) | rt.notation-cond-energies-geometry |  |
| `E_{C}` | upper characteristic energy E_ci + E_pe^2/E_ce | keV | `E_{C}` (rt.notation-cond-energies-geometry) | rt.notation-cond-energies-geometry |  |
| `\tilde E_{C}` | effective characteristic energy E_ci + (E~_pe)^2/E_ce at which the square of the effective refraction index becomes positive | keV | `\tilde E_{C}` (rt.notation-cond-energies-geometry) | rt.notation-cond-energies-geometry |  |
| `\tilde E_{pe}` | angle-dependent effective plasma energy E_pe sqrt(3 - 2 cos theta_k) | keV | `\tilde E_{pe}` (rt.notation-cond-energies-geometry) | rt.notation-cond-energies-geometry |  |
| `E_{L}` | line energy of the plasma-edge feature in the emissivity fit | keV | `E_{L}` (rt.notation-cond-energies-geometry) | rt.notation-cond-energies-geometry |  |
| `\theta_{kB}^{(i)}` | angle between incident wave vector and magnetic field, cos = sin theta_Bn sin theta cos phi_k - cos theta_Bn cos theta, with theta the polar angle of the folded incidence parameter (Pot12 theta_k) | rad | `\alpha_{i}` (rt.notation-cond-energies-geometry) | rt.notation-cond-energies-geometry |  |
| `\theta_{kB}^{(r)}` | angle between reflected (outgoing) wave vector and magnetic field, cos = sin theta_Bn sin theta cos phi_k + cos theta_Bn cos theta | rad | `\alpha_{r}` (rt.notation-cond-energies-geometry) | rt.notation-cond-energies-geometry |  |
| `\theta_{kB}^{(min)}` | min(theta_kB^(r), theta_kB^(i)); Pot12 writes the bare alpha for this angle | rad | `α` (rt.notation-cond-energies-geometry) | rt.notation-cond-energies-geometry |  |
| `J` | mean dimensionless emissivity J = 1 - R = (J_1 + J_2)/2 | 1 | `J` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `J_{l}` | dimensionless emissivity of linear polarisation l = 1 (perpendicular), 2 (parallel to the plane of incidence), J_l = 1 - R_l; orientation fixed by the explicit vectors of Pot12 App. B, which contradict the paper prose at lines 215-216 | 1 | `J_{j}` (rt.notation-cond-emissivity-fit); text/other: `J_1` | rt.notation-cond-emissivity-fit |  |
| `R_{mean}` | mean reflectivity (R_1 + R_2)/2 (paper: R) | 1 | `R` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `R_{l}` | effective reflectivity of linear polarisation j | 1 | `R_{j}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `r_{l'l}` | complex reflected-field amplitude into polarisation m for incident polarisation j | 1 | `r_{mj}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `\phi_{l'l}` | phase of r_mj | rad | `\phi_{mj}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `R_{jj'}` | reflection coefficient from incident mode M' into reflected mode M (M, M' = X, O), \|r'_MM'\|^2 | 1 | `R_{MM'}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit | Mode index M (Pot12) -> j. |
| `J_X` | normal-mode emissivity 1 - R_X (X-mode), with R_X = R_XX + R_XO | 1 | `J_{X}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `J_O` | normal-mode emissivity 1 - R_O (O-mode), with R_O = R_OO + R_OX | 1 | `J_{O}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `I_{E,l}` | specific intensity per unit energy emitted in linear polarisation l (l = 1 perpendicular, 2 parallel to the plane of incidence); I_E = (1.602176634e-9 erg/keV) I_nu/h for per-keV units | erg s^-1 cm^-2 keV^-1 sr^-1 | `I_{E,j}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit | l = 1 (in plane of incidence), 2 (perpendicular): linear basis, not X/O modes. |
| `I_E^{j}` | specific intensity per unit energy of normal mode M = X, O in the atmosphere; I_E = (1.602176634e-9 erg/keV) I_nu/h for per-keV units | erg s^-1 cm^-2 keV^-1 sr^-1 | `I_E^{M}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit | Per-energy intensity: I_E = I_nu/h. |
| `B_E` | Planck function per unit energy, B_E = (1.602176634e-9 erg/keV) B_nu/h for per-keV units | erg s^-1 cm^-2 keV^-1 sr^-1 | `B_E` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `P_{lin}` | degree of linear polarisation of the emitted radiation | 1 | `P_{lin}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `J_A` | fit emissivity in Region I (low energies) | 1 | `J_A` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `J_B` | fit emissivity for E_ci < E <~ E~_C (power-law interpolation) | 1 | `J_B` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `J_C` | fit factor for the steep rise at E > E~_C | 1 | `J_C` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `J_0` | approximate zero-angle emissivity 1 - (R_- + R_+)/2 of eq. 10, neglecting finite electron relaxation | 1 | `J_0` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `A_{fit}` | fit function A(E) of eq. 16 | 1 | `A(E)` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `L_{line}` | fit function describing the line near E_pe (paper: L) | 1 | `L_{line}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `X_L` | dimensionless line-profile variable (E - E_L)/(2 E_pe W_L)/(1 - cos theta_k) (paper: X) | 1 | `X` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `W_L` | line width parameter | 1 | `W_L` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `R_{L}` | line-strength function in the fit of linear polarisation l = 1 | 1 | `R_{L}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `J_{A1}` | fit term J_A of the linear-polarisation-1 emissivity J_1 = J_{l=1} | 1 | `J_{A1}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `J_{B1}` | fit term J_B of the linear-polarisation-1 emissivity J_1 = J_{l=1} | 1 | `J_{B1}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `a_1` | linear-polarisation-1 fit coefficient a_1 of eq. 25 | 1 | `a_1` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `A_1` | linear-polarisation-1 fit coefficient A_1 of eq. 25 | 1 | `A_1` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `p` | power-law exponent of J_B (fit of the mean emissivity) | 1 | `p` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `p_1` | power-law exponent of J_B1 (fit of linear polarisation 1) | 1 | `p_1` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `n_0` | refractive-index parameter (1 + E_pe^2/(2 E_ce E_ci))^(1/2) of eq. 19 | 1 | `n_0` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `\tilde n` | effective refraction index with n~^2 = 1 - E~_pe^2/(E_ce (E - E_ci)) | 1 | `\tilde n` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `n_\pm^{(0)}` | refraction indices for theta_Bn = theta = 0 (eq. 10), zero-angle approximation neglecting finite relaxation | 1 | `n_\pm^{(0)}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `R_\pm^{(0)}` | reflectivities for theta_Bn = theta = 0 (eq. 10), zero-angle approximation neglecting finite relaxation | 1 | `R_\pm^{(0)}` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `f_E` | smooth weight function E/(E + E~_C/2) of eqs. B.16-B.19 | 1 | `f_E` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `f_L` | logistic weight function of eqs. B.16-B.19 | 1 | `f_L` (rt.notation-cond-emissivity-fit) | rt.notation-cond-emissivity-fit |  |
| `\rho_{s,cond}` | zero-pressure density of the condensed surface (paper: rho_s) | g cm^-3 | `ρ_s` (rt.notation-cond-matter-dielectric); text/other: `rho_s` | rt.notation-cond-matter-dielectric |  |
| `\eta_{\rho}` | dimensionless factor eta in the surface density rho_s: uniform-electron (Wigner-Seitz) model, unity for the ion-sphere model | 1 | `η` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `T_{crit}` | critical temperature below which condensation occurs | K | `T_{crit}` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `Q_{s}` | cohesive energy of the condensed phase | eV | `Q_{s}` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `u_{e}` | (E_Be/E)^2, squared ratio of electron cyclotron energy to photon energy | 1 | `u_{e}` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `u_{i}` | (E_Bi/E)^2, squared ratio of ion cyclotron energy to photon energy | 1 | `u_{i}` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `v_{e}` | (E_pe/E)^2, squared ratio of electron plasma energy to photon energy | 1 | `v_{e}` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `\gamma_{ei}^{(tr)}` | dimensionless collisional damping rate hbar nu/E for motion transverse to B | 1 | `\gamma_{ei}^{(tr)}` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `\gamma_{ei}^{(l)}` | dimensionless collisional damping rate hbar nu/E for motion longitudinal to B | 1 | `\gamma_{ei}^{(l)}` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `\nu_{ei}^{(tr)}` | electron-ion collision frequency transverse to B | s^-1 | `\nu_{ei}^{(tr)}` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `\nu_{ei}^{(l)}` | electron-ion collision frequency longitudinal to B | s^-1 | `\nu_{ei}^{(l)}` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `\Gamma_{coul}` | Coulomb coupling parameter Gamma = (Ze)^2/(a_ion k_B T) of the ions (van Adelsberg+2005 sect. 2.1) | 1 | `Γ` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric |  |
| `a_{ion}` | ion-sphere radius a_i, defined by (4 pi a_i^3/3)^(-1) = n_i (ion number density) | cm | `a_i` (rt.notation-cond-matter-dielectric) | rt.notation-cond-matter-dielectric | Added in r1: the radius used by Gamma_coul. Not the BSk coefficients a_i^{BSk} or the magnetic length a_m. |

### Radiative transfer

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `\theta_{Bn}` | angle between the magnetic field B and the local surface normal n | rad | `\theta_B^{(n)}` (rt.notation-cond-energies-geometry); `\Theta_B` (obs.notation-rxj1856-spectral); `θ_Bn` (rt.notation-spw09); text/other: `Theta` | rt.notation-spw09 | Ho07 tables quote degrees. |
| `\theta` | polar angle between the photon wave vector k and the local outward surface normal n; mu = cos theta (Pot12 theta_k is a folded incidence parameter: it equals theta for outgoing rays and pi - theta for incident rays) | rad | `\theta_k` (rt.notation-cond-energies-geometry); `θ` (rt.notation-spw09); text/other: `theta^(i)` | rt.notation-spw09 |  |
| `\varphi_k` | azimuthal angle of the photon wave vector about the surface normal, measured from the plane containing B and the normal | rad | `\varphi_{k}` (rt.notation-cond-energies-geometry); `φ` (rt.notation-spw09); text/other: `phi` | rt.notation-spw09 |  |
| `\mu` | cosine of the angle between a ray and the outward surface normal (positive); ingoing rays are written with -μ | 1 | `μ` (rt.notation-mcphac) | rt.notation-mcphac | Cosine of the ray angle to the outward normal, mu = cos theta; McPHAC half-range (> 0, ingoing written -mu), SPW09/VAL06 signed. Never cos psi, a chemical potential or a reduced mass. (r1 fix: adopt signed mu = cos theta in [-1, 1] project-wide; McPHAC and VAL06 both use a positive half-range argument with explicit +-mu.) |
| `I_\nu(m,\mu)` | specific intensity of outgoing radiation at column depth y, frequency ν, direction μ | erg s^-1 cm^-2 Hz^-1 sr^-1 | `I_ν(y,μ)` (rt.notation-mcphac) | rt.notation-mcphac | (r1 fix blocked: an evaluation of the existing I_nu; rt.notation-mcphac is blessed and needs the call-form row, so the row is retained; owner decision needed. No edit.) |
| `I_\nu(m,-\mu)` | specific intensity of ingoing radiation (direction -μ) | erg s^-1 cm^-2 Hz^-1 sr^-1 | `I_ν(y,-μ)` (rt.notation-mcphac) | rt.notation-mcphac | (r1 fix blocked: folding into the signed-direction I_nu(m,mu) would delete a row of the blessed rt.notation-mcphac and change its symbol set; owner decision needed. No edit.) |
| `P_\nu(m,\mu)` | Feautrier variable: half the sum of outgoing and ingoing specific intensity, (I(y,μ)+I(y,-μ))/2 | erg s^-1 cm^-2 Hz^-1 sr^-1 | `P_ν(y,μ)` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `m` | column density (mass per unit area above the point), dm = -rho dz | g cm^-2 | `y_col` (rt.notation-mcphac); `m` (rt.notation-spw09); text/other: `y` | rt.notation-spw09 | McPHAC and pipeline entries write y; y is reserved for CIE chromaticity. |
| `P_g` | gas pressure (electrons plus ions) | dyn cm^-2 | `P(y)` (rt.notation-mcphac); `P_g` (rt.notation-spw09); `P_gas` (rt.notation-val); text/other: `P` | rt.notation-spw09 | Bare P is the total (EOS/TOV) pressure; equal to P_g in the atmosphere when radiation pressure is neglected. |
| `\tau_\nu` | effective normal optical depth, dτ_ν = k_ν dy | 1 | `τ_ν` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `\varpi_\nu` | dimensionless scattering albedo σ_T,r/(κ_ν+σ_T,r) = σ_T,r/k_ν (not mass density) | 1 | `ρ_ν` (rt.notation-mcphac) | rt.notation-mcphac | Albedo; not a density. |
| `B_\nu(T)` | blackbody (Planck) intensity at temperature T | erg s^-1 cm^-2 Hz^-1 sr^-1 | `B_ν(T)` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `S_\nu(m,\mu)` | source function | erg s^-1 cm^-2 Hz^-1 sr^-1 | `S_ν(y,μ)` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `d\sigma_{T,r}/d\mu'(\mu,\mu')` | differential amount of scattering from directions μ′ into directions μ | cm^2 g^-1 | `dσ_T,r/dμ′(μ,μ′)` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `\mu_q` | Gauss-Legendre angle nodes on (0,1); M of them | 1 | `μ_j` (rt.notation-mcphac) | rt.notation-mcphac | Gauss-Legendre node index q (McPHAC writes j, which is the mode index here). |
| `a_q` | Gauss-Legendre weights on (0,1); M of them | 1 | `a_j` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `J_\nu` | mean specific intensity, integral over 0..1 of P_nu dmu; the moment-equation J_nu equals the transfer-solution mean only as the temperature correction vanishes (McPHAC) | erg s^-1 cm^-2 Hz^-1 sr^-1 | `J_ν` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `f_\nu` | Eddington factor (second moment of P_ν over J_ν) | 1 | `f_ν` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `h_\nu` | surface Eddington factor (first moment of P_ν over J_ν) | 1 | `h_ν` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `k_R` | Rosseland mean total absorption opacity | cm^2 g^-1 | `k_R` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `T^{(0)}` | present best-guess temperature profile T_0(y) of the iteration, T = T_0 + Delta T | K | `T0(y)` (rt.notation-mcphac); text/other: `T0`, `T_0(y)` | rt.notation-mcphac |  |
| `\Delta T` | temperature correction Delta T(y), T(y) = T_0(y) + Delta T(y) | K | `ΔT` (rt.notation-mcphac) | rt.notation-mcphac |  |
| `I_\nu^j` | specific intensity in normal mode j (j = 1 X, 2 O) as a function of depth and direction; +mu outgoing, -mu ingoing | erg s^-1 cm^-2 Hz^-1 sr^-1 | `I_ν^i(μ,φ)` (rt.notation-spw09); `I_ν^j(τ,±μ)` (rt.notation-val); text/other: `I^(j)_nu`, `I_O`, `I_X` | rt.notation-val |  |
| `z_{vert}` | vertical distance z in the atmosphere (SPW09) | cm | `z_vert` (rt.notation-spw09) | rt.notation-spw09 |  |
| `g_{rad}` | radiative acceleration | cm s^-2 | `g_rad` (rt.notation-spw09) | rt.notation-spw09 |  |
| `d\kappa^s_{ji}/d\Omega` | differential scattering opacity from mode j direction (mu',phi') into mode i direction (mu,phi); under SPW09 definitions dkappa^s_{ji}/dOmega = sigma_nu^{ij}/(2 pi), not sigma_nu^{ij} itself | cm^2 g^-1 sr^-1 | `σ_ν^ij(μ,φ;μ′,φ′)` (rt.notation-spw09) | rt.notation-spw09 |  |
| `\tau_\nu^j` | optical depth in mode i, dτ_ν^i = (k_ν^i+σ_ν^i) dm | 1 | `τ_ν^i` (rt.notation-spw09) | rt.notation-spw09 |  |
| `S_\nu^j` | source function of normal mode j | erg s^-1 cm^-2 Hz^-1 sr^-1 | `S_ν^i` (rt.notation-spw09); `S_ν^j(μ)` (rt.notation-val) | rt.notation-val |  |
| `\eta_\nu^j` | emissivity in mode i, per unit mass-column normalisation of eq. 10 | erg s^-1 g^-1 Hz^-1 sr^-1 | `η_ν^i` (rt.notation-spw09) | rt.notation-spw09 |  |
| `H_\nu` | Eddington flux at column density m | erg s^-1 cm^-2 Hz^-1 | `H_ν(m)` (rt.notation-spw09) | rt.notation-spw09 |  |
| `H_0` | σT_eff^4/(4π), the required Eddington flux integrated over frequency | erg s^-1 cm^-2 | `H_0` (rt.notation-spw09) | rt.notation-spw09 |  |
| `\varepsilon_H` | relative flux error 1 - H_0/∫H_ν dν | 1 | `ε_H(m)` (rt.notation-spw09) | rt.notation-spw09 |  |
| `\varepsilon_\Lambda` | energy-balance error as a function of depth | erg s^-1 g^-1 | `ε_Λ(m)` (rt.notation-spw09) | rt.notation-spw09 |  |
| `\Lambda^j_{\nu,diag}` | diagonal element of the Λ operator for mode i, a function of the mean optical depth | 1 | `Λ^i_ν,diag` (rt.notation-spw09) | rt.notation-spw09 |  |
| `\varpi_\nu^j` | scattering albedo of mode i from averaged opacities, σ̄/(k̄+σ̄) | 1 | `α_ν^i` (rt.notation-spw09) | rt.notation-spw09 |  |
| `\Delta T_\Lambda` | temperature correction from the integral Λ-iteration | K | `ΔT_Λ` (rt.notation-spw09) | rt.notation-spw09 |  |
| `\tau_T` | Thomson optical depth, dτ = -ρ κ_T dz | 1 | `τ` (rt.notation-val) | rt.notation-val |  |
| `J_\nu^j` | mean intensity of mode j, (1/2)∫_0^1 dμ [I(μ)+I(-μ)] | erg s^-1 cm^-2 Hz^-1 sr^-1 | `J_ν^j` (rt.notation-val) | rt.notation-val |  |
| `u_\nu^j` | specific energy density of mode j, (2π/c)∫_-1^1 dμ′ I_ν^j(μ′) | erg cm^-3 Hz^-1 | `u_ν^j` (rt.notation-val) | rt.notation-val |  |
| `B_\lambda` | blackbody (Planck) intensity per unit wavelength, B_lambda = B_nu c/lambda^2 | erg s^-1 cm^-3 sr^-1 (per cm; per nm with 1 nm = 1e-7 cm) | `B_lambda` (rt.planck-function-blambda; Marr & Wilkin 2011) | rt.planck-function-blambda | (M1 2026-10-23) Companion of `B_\nu(T)`; not B_nu with lambda = c/nu substituted (the Jacobian c/lambda^2 is part of it). |

### Surface model and observer-frame quantities

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `T_s` | local effective surface temperature in the local reference frame at the stellar surface | K | `T_s` (eos.notation-structure, surface.notation-temperatures); `T_eff,atm` (rt.notation-mcphac) | surface.notation-temperatures | (r1 fix: McPHAC T_eff,atm is a contextual alias of the local T_s.) |
| `T_s^\infty` | redshifted (apparent) surface temperature as seen by a distant observer | K | `T_s^inf` (eos.notation-structure) | eos.notation-structure | Local; global analogue is T_eff^inf. (r1 fix applied: potekhin_yakovlev_2001_cooling.txt L227-233 "Ts∞ = Ts 1 − rg /R" (sqrt sign lost by the text layer); recorded in notes of truth/entries/eos.toml:43; no per-symbol source field.) |
| `L_{ph}` | photon luminosity in the local reference frame at the stellar surface | erg s^-1 | `L` (eos.notation-structure); `L_ph` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `L_{ph}^\infty` | apparent photon luminosity for a distant observer, L_ph^inf = (1 - r_g/R) L_ph | erg s^-1 | `L^inf` (eos.notation-structure); `L_ph^inf` (gr.notation-potekhin); `L_ph^∞` (surface.notation-observer-frame) | surface.notation-observer-frame |  |
| `R^\infty` | apparent radius of the star for a distant observer, R^inf = R/sqrt(1 - r_g/R) = R(1 + z_g) | cm | `R^inf` (eos.notation-structure); `R_inf` (gr.notation-potekhin); `R^∞` (surface.notation-observer-frame); text/other: `R_infty` | surface.notation-observer-frame |  |
| `T_{eff}^\infty` | redshifted (apparent) effective temperature seen by a distant observer, T_eff^inf = T_eff sqrt(1 - r_g/R) = T_eff/(1 + z_g) | K | `T_eff^inf` (gr.notation-potekhin); `T^\infty` (obs.notation-rxj1856-spectral); `T_eff^∞` (surface.notation-temperatures) | surface.notation-temperatures | (disputed: r1 rejected merging the global flux-defined T_eff^inf (PPP15:2098-2111) with the fitted emission-region temperature of Ho07 (ho_2007_rxj1856_thin_atm:85-88,345-347); the intended split is T_eff^inf (global) and a separately labelled T_em^inf (fitted).) |
| `z_g` | gravitational redshift, 1 + z_g = (1 - r_g/R)^(-1/2) = (1 - x_g)^(-1/2) | 1 | `z_g^{fit}` (obs.notation-rxj1856-spectral); `z_g` (surface.notation-observer-frame); text/other: `z` | surface.notation-observer-frame |  |
| `\theta_B` | magnetic colatitude: polar angle on the stellar surface measured from the magnetic axis (written beta in GH83 eq. 15 and theta in PPP15 App. B) | rad | `θ_B` (surface.notation-geometry); text/other: `beta (GH83)`, `theta (PPP15)` | surface.notation-geometry | Magnetic colatitude only (CLAUDE.md usage). The k-B angle is theta_kB; the field-normal angle is theta_Bn. |
| `\Phi_B` | angle between the magnetic field lines and the local horizontal (surface tangent plane) (GH83 Fig. 1, Phi) | rad | `Φ` (surface.notation-geometry); text/other: `Φ_B` | surface.notation-geometry | (r1 fix applied: units rad; source spelling Phi.) |
| `\Psi_{pole}` | geographic colatitude (angle from the rotation axis) of a magnetic pole (GH83 Fig. 3, Psi) | rad | `Ψ` (surface.notation-geometry); text/other: `Ψ_pole` | surface.notation-geometry |  |
| `\theta_{los}` | colatitude of the line of sight to the observer measured from the rotation axis (written theta in GH83 Fig. 3) | rad | `θ` (surface.notation-geometry); text/other: `θ_los` | surface.notation-geometry |  |
| `\phi_{rot}` | pulsar rotational phase, phi = Omega t (GH83 eq. 21 and Fig. 3 caption) | rad | `φ` (surface.notation-geometry); text/other: `φ_rot` | surface.notation-geometry |  |
| `\psi_{pole}` | angle between the normal to a magnetic pole and the line of sight to the observer (GH83 Sect. III a; written alpha) | rad | `α` (surface.notation-geometry); text/other: `α_pole` | surface.notation-geometry | GH83 alpha: angle between the pole normal and the line of sight = psi of Bel02 evaluated at the pole. |
| `T_{eff}` | overall effective temperature of the star in the local frame, 4 pi sigma_SB R^2 T_eff^4 = L_ph = integral of F_ph dSigma (PPP15 eq. 39) | K | `T_eff` (surface.notation-temperatures) | surface.notation-temperatures |  |
| `T_{pole}` | surface temperature at the magnetic pole (T_pole in GH83 eq. 14; T_p in PPP15 App. B) | K | `T_pole` (surface.notation-temperatures) | surface.notation-temperatures |  |
| `T_{eq}` | surface temperature at the magnetic equator (PPP15 App. B, T_eq) | K | `T_eq` (surface.notation-temperatures) | surface.notation-temperatures |  |
| `T_b` | temperature at the bottom of the heat-blanketing envelope (crust/envelope boundary, density rho_b) (PPP15 Sect. 5.1) | K | `T_b` (surface.notation-temperatures) | surface.notation-temperatures | (r1 fix applied: locator of surface.notation-temperatures now cites Sect. 5.1 (potekhin_pons_page_2015_thermal.txt L1942-1943 "temperature Tb at the crust/envelope boundary") and App. B L2674; truth/entries/surface.toml:50.) |
| `T_*` | temperature from the Appendix-B Ts fits evaluated without the neutrino-limit correction (B.4), used to define the flux F_b = sigma_SB T_*^4 (PPP15 App. B item 4) | K | `T_*` (surface.notation-temperatures) | surface.notation-temperatures |  |
| `B_p` | magnetic field strength at the magnetic pole (PPP15 App. B) | G | `B_p` (surface.notation-envelope-fit); text/other: `B_pole` | surface.notation-envelope-fit | Polar field. The spin-down B_dip is a different (equatorial-convention) estimate. (r1 fix: polar field only; no equatorial convention is asserted for the separate spin-down estimate B_dip.) |
| `B_{p,12}` | polar magnetic field B_p/(10^12 G) in Appendix B (PPP15 App. B: B12 = Bp/1012 G) | 1 | `B_12` (surface.notation-envelope-fit); text/other: `B12 (PPP15 App. B)` | surface.notation-envelope-fit |  |
| `T_9` | T_b / 1e9 K | 1 | `T_9` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `g_{14}` | surface gravity in units of 1e14 cm s^-2 (used in PPP15 eqs. 37, B.1, B.3; not defined explicitly in PPP15 itself) | 1 | `g_14` (surface.notation-envelope-fit) | surface.notation-envelope-fit | (r1 fix applied: missing-file citation replaced by potekhin_yakovlev_2001_cooling.txt L539-540 "g14 ≡ g/1014 cm s−2"; truth/entries/surface.toml:86.) |
| `T_{pole}^{(0)}` | polar surface temperature from fit B.1, neglecting crustal neutrino emission (T_p^(0) in the paper) | K | `T_pole^(0)` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `T_{pole}^{(max)}` | limiting polar surface temperature where the Tp(Tb) relation levels off due to crustal neutrino emission, eq. B.3 | K | `T_pole^(max)` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `T_0^{B.2}` | dimensionless fit variable T_0 defined in eq. B.2 (not the GH83 reservoir temperature and not a dimensional temperature scale) | 1 | `T_0` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `T_1^{B.2}` | dimensionless fit variable T_1 defined in eq. B.2 | 1 | `T_1` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `a^{B.2}` | exponent a = 0.337/(1 + 0.02 sqrt(B_{p,12})) in eq. B.2, with B_{p,12} the polar field in 1e12 G | 1 | `a` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `a_1^{B.6}` | coefficient a_1 = a_2 T_9^{1/2}/3 in eq. B.6 | 1 | `a_1^{B.6}` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `a_2^{B.6}` | coefficient a_2 = 10 B_{p,12}/(T_9^{1/2} + 0.1 B_{p,12} T_9^{-1/4}) in eq. B.6, with B_{p,12} the polar field in 1e12 G | 1 | `a_2^{B.6}` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `\rho_b` | mass density at the bottom boundary of the heat-blanketing envelope (1e10 g cm^-3 for the Appendix-B fit) | g cm^-3 | `ρ_b` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `F_b` | heat flux density through the envelope bottom boundary rho = rho_b | erg cm^-2 s^-1 | `F_b` (surface.notation-envelope-fit) | surface.notation-envelope-fit |  |
| `F_{ph}` | local radial photon energy flux density at the stellar surface | erg cm^-2 s^-1 | `F_rad` (rt.notation-val); `F_ph` (surface.notation-envelope-fit) | surface.notation-envelope-fit | Equals F_rad of the atmosphere at that patch. (r1 fix: VAL06 F_rad at the atmosphere surface is identified with F_ph.) |
| `K_\parallel` | thermal conductivity along the magnetic field B (GH83 eq. 6) | erg cm^-1 s^-1 K^-1 | `K_∥` (surface.notation-conduction) | surface.notation-conduction |  |
| `K_\perp` | thermal conductivity perpendicular to the magnetic field B (GH83 eq. 6); the anisotropic model assumes K_perp << K_par, the source also considers equality | erg cm^-1 s^-1 K^-1 | `K_⊥` (surface.notation-conduction) | surface.notation-conduction |  |
| `K_{ab}^{cond}` | thermal conductivity tensor K_ab, spatial indices a, b in {x, z} in the horizontal/vertical frame (GH83 eqs. 7-8); not the mode indices i, j | erg cm^-1 s^-1 K^-1 | `K_ij` (surface.notation-conduction) | surface.notation-conduction | (r1 fix: renamed from K_{ij} so that spatial conductivity indices do not collide with the mode indices i, j.) |
| `q_z` | upward heat-flux component at the surface, q_z(surface) = sigma T^4 (GH83 eqs. 10-11) | erg cm^-2 s^-1 | `q_z` (surface.notation-conduction) | surface.notation-conduction |  |
| `L_{slab}` | depth of the isothermal reservoir below the surface (written L in GH83 Fig. 1) | cm | `L` (surface.notation-conduction); text/other: `L_slab` | surface.notation-conduction |  |
| `R_{cap}` | radius of the polar cap, R_cap = R^(3/2) [2 pi/(c P_spin)]^(1/2) << R (GH83 eq. 1, where the period is written p) | cm | `R_cap` (surface.notation-conduction) | surface.notation-conduction |  |
| `M_{mod}` | degree of modulation: (peak intensity - minimum intensity)/minimum intensity (GH83 eq. 23) | 1 | `M` (surface.notation-conduction); text/other: `M_mod` | surface.notation-conduction |  |
| `T_{cap}` | temperature of the heated polar cap (GH83 eq. 2-3) | K | `T_cap` (surface.notation-conduction) | surface.notation-conduction |  |
| `T_{bg}` | background temperature of the stellar surface outside the cap (written T in GH83 eqs. 2-3) | K | `T` (surface.notation-conduction); text/other: `T_bg` | surface.notation-conduction |  |
| `T_{res}` | temperature of the isothermal heat reservoir below the surface slab (written T_0 in GH83 Fig. 1) | K | `T_0` (surface.notation-conduction); text/other: `T_res` | surface.notation-conduction |  |
| `T_{fit}` | best-fit blackbody temperature of a synthesised spectrum (GH83 Sect. III c; the source describes it in prose, T_fit is project notation) | K | `T_fit` (surface.notation-conduction) | surface.notation-conduction |  |
| `a_{dip}` | GR dipole parameter a(x), x_g = r_g/R, in tan theta_Bn = a tan theta_B and B = B_p sqrt(cos^2 theta_B + a^2 sin^2 theta_B); a -> 1/2 as x_g -> 0 | 1 | `a` (surface.py01-gr-dipole-a-of-x, surface.py01-gr-dipole-magnitude, surface.py01-gr-dipole-field-angle) | surface.py01-gr-dipole-a-of-x | (M1 2026-10-23) Not the other a coefficients (`a_m`, `a_1^{B.6}`, `a_{key}`, `a_{spin}`). |
| `\theta_{surf}` | colatitude of a surface point from the rotation axis, in the non-rotating frame with z along the rotation axis | rad | (surface.field-colatitude-from-geometry) | surface.field-colatitude-from-geometry | (M1 2026-10-23) Project convention. Not `\theta_{Sch}` (Bel02 orbit-plane polar coordinate). |
| `\phi_{surf}` | azimuth of a surface point about the rotation axis in the non-rotating frame; a material point has phi_surf = phi_0 + phi_rot | rad | (surface.field-colatitude-from-geometry) | surface.field-colatitude-from-geometry | (M1 2026-10-23) Project convention. |
| `\mu_{mag}` | magnetic dipole moment of the star (WS83 mu) | G cm^3 | `mu` (surface.ws83-gr-dipole-components; WS83 eqs. 12a-b) | surface.ws83-gr-dipole-components | (M1 2026-10-23) Not the direction cosine `\mu` or the Bohr magneton `\mu_B`. |

### Observational (RX J1856.5-3754)

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `R_{em}^\infty` | apparent (radiation) radius at infinity, R^inf = Rem (1+z_g) | km | `R^\infty` (obs.notation-rxj1856-spectral) | obs.notation-rxj1856-spectral |  |
| `R_{em}` | physical size of the emission region | km | `R_{em}` (obs.notation-rxj1856-spectral) | obs.notation-rxj1856-spectral |  |
| `N_H` | interstellar hydrogen column density along the line of sight to the source | cm^-2 | `N_H^{fit}` (obs.notation-rxj1856-spectral); `N_H` (surface.notation-conduction) | obs.notation-rxj1856-spectral |  |
| `y_H` | hydrogen column depth of the thin atmosphere (mass per area) | g cm^-2 | `y_H` (obs.notation-rxj1856-spectral) | obs.notation-rxj1856-spectral |  |
| `F_\lambda` | flux density per unit wavelength | erg cm^-2 s^-1 A^-1 | `F_\lambda` (obs.notation-rxj1856-photometry) | obs.notation-rxj1856-photometry |  |
| `\alpha_{opt}` | (paper symbol alpha) optical/UV power-law index, F_lambda ~ lambda^-alpha (Rayleigh-Jeans: 4) | 1 | `\alpha_{opt}` (obs.notation-rxj1856-photometry) | obs.notation-rxj1856-photometry |  |
| `A_V` | visual extinction | mag | `A_V` (obs.notation-rxj1856-photometry) | obs.notation-rxj1856-photometry |  |
| `R^\infty_{bb}/D` | X-ray blackbody normalisation (apparent radius per distance) | km kpc^-1 | `R/d` (obs.notation-rxj1856-photometry) | obs.notation-rxj1856-photometry |  |
| `\mathrm{STMAG}` | magnitude in the STMAG system | mag | `\mathrm{STMAG}` (obs.notation-rxj1856-photometry) | obs.notation-rxj1856-photometry |  |
| `P_{spin}` | spin (rotation) period of the neutron star | s | `P_{spin}` (obs.notation-rxj1856-timing); `P_spin` (surface.notation-conduction); text/other: `P (vK&K08)`, `P_{1856}`, `p (GH83)` | obs.notation-rxj1856-timing |  |
| `\dot P` | period derivative | s s^-1 | `\dot P` (obs.notation-rxj1856-timing) | obs.notation-rxj1856-timing |  |
| `\nu_{spin}` | spin frequency (paper symbol nu) | Hz | `\nu_{spin}` (obs.notation-rxj1856-timing) | obs.notation-rxj1856-timing |  |
| `\dot\nu_{spin}` | spin frequency derivative (paper symbol nu-dot) | Hz s^-1 | `\dot\nu_{spin}` (obs.notation-rxj1856-timing) | obs.notation-rxj1856-timing |  |
| `B_{dip}` | dipole magnetic field inferred from spin-down | G | `B_{dip}` (obs.notation-rxj1856-timing) | obs.notation-rxj1856-timing |  |
| `\tau_{char}` | characteristic age P/(2 Pdot) | Myr | `\tau_{char}` (obs.notation-rxj1856-timing) | obs.notation-rxj1856-timing |  |

### Astrometry and Gaia photometry (M1 2026-10-23)

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `\alpha_{RA}` | right ascension of a source; the reference frame and epoch are attributes of each datum (Gaia DR3: ICRS, Gaia-CRF3, epoch J2016.0; Walter+2010: not stated) | rad | `alpha` (obs.notation-gaia-sky, obs.rxj1856-position-walter2010); text/other: `ra` | obs.notation-gaia-sky | (M1 2026-10-23) Catalogue units deg (Gaia) or h m s (Walter+2010); converted to rad. Bare alpha is the emission angle. Meaning generalised from the Gaia-only text so that the Walter+2010 datum shares the row (entry meaning edit pending). |
| `\delta_{Dec}` | declination of a source; the reference frame and epoch are attributes of each datum (Gaia DR3: ICRS, Gaia-CRF3, epoch J2016.0; Walter+2010: not stated) | rad | `delta` (obs.notation-gaia-sky, obs.rxj1856-position-walter2010); text/other: `dec` | obs.notation-gaia-sky | (M1 2026-10-23) Catalogue units deg or deg arcmin arcsec. Bare delta is the Pech83 alias of the emission angle; `\delta_1` is the DA09 geodesic switch. Meaning edit pending as for `\alpha_{RA}`. |
| `\varpi_{plx}` | trigonometric parallax of a source at the reference epoch J2016.0 (Gaia astrometric parameter) | rad | `varpi` (obs.notation-gaia-sky, obs.gaia-perspective-acceleration) | obs.notation-gaia-sky | (M1 2026-10-23) Catalogue unit mas. Subscripted to keep clear of the albedos `\varpi_\nu`, `\varpi_\nu^j`. |
| `\mu_{\alpha*}` | proper motion in right ascension times cos(delta), at J2016.0 (catalogue field pmra) | rad s^-1 | `mu_{alpha*} = mu_alpha cos(delta)` (obs.notation-gaia-sky) | obs.notation-gaia-sky | (M1 2026-10-23) Catalogue unit mas yr^-1. Not the ray cosine `\mu`. |
| `\mu_\delta` | proper motion in declination at J2016.0 (catalogue field pmdec) | rad s^-1 | `mu_delta` (obs.notation-gaia-sky) | obs.notation-gaia-sky | (M1 2026-10-23) Catalogue unit mas yr^-1. Not the ray cosine `\mu`. |
| `v_r` | radial velocity of a source relative to the solar system barycentre (sixth kinematic parameter of the standard model) | cm s^-1 | `v_r` (obs.notation-gaia-sky) | obs.notation-gaia-sky | (M1 2026-10-23) Catalogue unit km s^-1. |
| `G_{Gaia}` | Gaia G-band magnitude | mag | `G` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | (M1 2026-10-23) Bare G is Newton's constant. |
| `G_{BP}` | Gaia blue-photometer band magnitude | mag | `G_BP` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | (M1 2026-10-23) |
| `G_{RP}` | Gaia red-photometer band magnitude | mag | `G_RP` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | (M1 2026-10-23) |
| `n_{p,Gaia}` | mean source flux in a Gaia band in electrons per second (e- s^-1), not normalised by telescope area; in the XP chapter n_p(lambda) is the photon flux distribution | s^-1 | `n_p` (obs.notation-gaia-photometry, obs.gaia-band-flux-integral) | obs.notation-gaia-photometry | (M1 2026-10-23) Bare `n_p` is the proton number density. |
| `S_{Gaia}` | Gaia system overall response function, including the scaling from photons to electrons (passband) | e- photon^-1 | `S(lambda)` (obs.notation-gaia-photometry, obs.gaia-band-flux-integral, obs.gaia-synthetic-mean-fluxes, obs.gaia-vegamag-definition) | obs.notation-gaia-photometry | (M1 2026-10-23) Bare S is retired. |
| `P_{pupil}` | telescope pupil area in the Gaia flux integral | cm^2 | `P` (obs.notation-gaia-photometry, obs.gaia-band-flux-integral) | obs.notation-gaia-photometry | (M1 2026-10-23) Printed SI unit m^2 (units_source). Bare P is pressure. |
| `ZP_{VEG}` | Gaia passband zero point in the VEGAMAG system, per band | mag | `ZP_{VEG}` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | (M1 2026-10-23) |
| `ZP_{AB}` | Gaia passband zero point in the AB system, per band | mag | `ZP_{AB}` (obs.notation-gaia-photometry) | obs.notation-gaia-photometry | (M1 2026-10-23) |
| `\langle F_\lambda\rangle` | synthetic mean energy flux per wavelength unit through a passband | erg s^-1 cm^-2 cm^-1 | `<f_lambda>` (obs.notation-gaia-photometry, obs.gaia-synthetic-mean-fluxes, obs.gaia-vegamag-definition) | obs.notation-gaia-photometry | (M1 2026-10-23) Renamed from the proposed `\langle f_\lambda\rangle`: flux density is `F_\lambda` in this dictionary and lower-case `f_\nu` is the Eddington factor (sky r2 decision, applied in the entry). Printed SI unit in units_source. |
| `\langle F_\nu\rangle` | synthetic mean energy flux per frequency unit through a passband | erg s^-1 cm^-2 Hz^-1 | `<f_nu>` (obs.notation-gaia-photometry, obs.gaia-synthetic-mean-fluxes) | obs.notation-gaia-photometry | (M1 2026-10-23) Renamed from `\langle f_\nu\rangle` for the same reason; pairs with the reserved `F_\nu` (sky r2 decision, applied in the entry). Printed SI unit in units_source. |
| `t_{ep}` | reference epoch of the astrometric parameters (arbitrary, chosen near the mid-time of the mission); Gaia DR3: J2016.0 | s | `t_ep` (obs.notation-gaia-propagation, obs.agis-standard-model-direction) | obs.notation-gaia-propagation | (M1 2026-10-23) |
| `t_{pr}` | time elapsed since the reference epoch over which the astrometric parameters are propagated (Julian years of 365.25 x 86400 s, TCB) | s | `t = T - T_0` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) Catalogue unit yr. Bare `t` is Schwarzschild coordinate time. |
| `\vec{u}_{bary}` | barycentric coordinate direction to the source, unit vector b/\|b\| in the ICRS | 1 | `u` (obs.notation-gaia-propagation, obs.bl-proper-motion-components, obs.bl-propagation-direction-nolt, obs.bl-propagation-direction-rigorous) | obs.notation-gaia-propagation | (M1 2026-10-23) Bare u is retired. |
| `\vec{r}_{tri}` | fixed unit vector coinciding with the barycentric direction at the reference epoch; third axis of the normal triad | 1 | `r` (obs.notation-gaia-propagation, obs.agis-normal-triad) | obs.notation-gaia-propagation | (M1 2026-10-23) Not the radial coordinate `r`. |
| `\vec{p}_{tri}` | unit vector in the direction of increasing right ascension (first axis of the normal triad) | 1 | `p` (obs.notation-gaia-propagation, obs.agis-normal-triad) | obs.notation-gaia-propagation | (M1 2026-10-23) Not the Pot12 exponent `p`. |
| `\vec{q}_{tri}` | unit vector in the direction of increasing declination (second axis of the normal triad) | 1 | `q` (obs.notation-gaia-propagation, obs.agis-normal-triad) | obs.notation-gaia-propagation | (M1 2026-10-23) Not the node index q or `q_{vac}`. |
| `\vec{\mu}_{pm}` | proper motion vector du/dt of the barycentric direction | rad s^-1 | `mu (vector)` (obs.notation-gaia-propagation, obs.bl-proper-motion-components) | obs.notation-gaia-propagation | (M1 2026-10-23) Catalogue unit mas yr^-1. |
| `\mu_{r}` | radial proper motion mu_r = v_r varpi / A (radial velocity in the unit of proper motion), sixth astrometric parameter | rad s^-1 | `mu_r` (obs.notation-gaia-propagation, obs.gaia-perspective-acceleration) | obs.notation-gaia-propagation | (M1 2026-10-23) Catalogue unit mas yr^-1. A is `\mathrm{au}`. |
| `\mu_{tot}` | total tangential proper motion, (mu_alpha*^2 + mu_delta^2)^(1/2) = v_t / b | rad s^-1 | `mu` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) Catalogue unit mas yr^-1. |
| `\vec{b}_{bary}` | barycentric position vector of the source (apparent unless stated) | cm | `b (vector)` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) Not the impact parameter `b`. |
| `b_{bary}` | barycentric distance of the source, varpi = A/b | cm | `b = \|b\|` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) Distinct from `D` (observer-star distance); the two differ by the observer's offset from the barycentre. |
| `\vec{v}_{bary}` | space velocity of the source relative to the solar system barycentre (apparent unless stated) | cm s^-1 | `v (vector)` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) |
| `v_{t}` | apparent tangential velocity of the source, v_t = mu b | cm s^-1 | `v_t` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) |
| `A_{vel}` | astronomical unit expressed in km yr s^-1: A_V = 4.740 470 446 | km yr s^-1 | `A_V` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) Bare `A_V` is the visual extinction. A conversion factor derived from `\mathrm{au}`. |
| `\tau_{A}` | light-travel time of one astronomical unit, tau_A = A/c = 499.004 784 s | s | `tau_A` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) |
| `\tau_{lt}` | light-travel time T_obs - T_em from the source to the barycentre; tau_0 = b_0/c at the reference epoch | s | `tau` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) Bare tau is taken by optical depths. |
| `f_{D}` | distance factor b_0/b of the epoch propagation | 1 | `f_D` (obs.notation-gaia-propagation, obs.bl-distance-factor-nolt, obs.bl-propagation-direction-nolt) | obs.notation-gaia-propagation | (M1 2026-10-23) |
| `f_{T}` | time factor of the epoch propagation with light-time effects | 1 | `f_T` (obs.notation-gaia-propagation, obs.bl-propagation-direction-rigorous) | obs.notation-gaia-propagation | (M1 2026-10-23) |
| `f_{V}` | velocity factor of the epoch propagation with light-time effects | 1 | `f_V` (obs.notation-gaia-propagation, obs.bl-propagation-pm-rigorous) | obs.notation-gaia-propagation | (M1 2026-10-23) |
| `\varepsilon_{T}` | small parameter f_T - 1 of the light-time expansion | 1 | `epsilon_T` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) |
| `\varepsilon_{V}` | small parameter f_V - 1 of the light-time expansion | 1 | `epsilon_V` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) |
| `\Delta\theta_{lt}` | angular-position error caused by neglecting light-time effects in the epoch propagation (first-order estimate) | rad | `Delta theta` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | (M1 2026-10-23) |

### Colour and display

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `L_w` | world luminance of pixel (x_pix,y_pix), scene luminance before scaling; source-unspecified units (absolute CGS/stilb is not asserted) | source-unspecified (luminance) | `L_w` (colour.notation-reinhard-luminance) | colour.notation-reinhard-luminance |  |
| `\bar{L}_w` | log-average luminance of the scene (eq. 1), same units as L_w | source-unspecified (luminance) | `\bar{L}_w` (colour.notation-reinhard-luminance) | colour.notation-reinhard-luminance |  |
| `L_{sc}(x_{pix},y_{pix})` | scaled luminance after key scaling, eq. 2 | 1 | `L(x,y)` (colour.notation-reinhard-luminance) | colour.notation-reinhard-luminance |  |
| `L_d` | display luminance output of the tone-mapping operator, eq. 3-4; eq. 4 permits values above 1 when L_sc exceeds L_white | 1 | `L_d` (colour.notation-reinhard-luminance) | colour.notation-reinhard-luminance |  |
| `L_{max}` | maximum scaled luminance in the scene, L_max = max L_sc (for use with eq. 4, alongside L_white) | 1 | `L_{max}` (colour.notation-reinhard-luminance) | colour.notation-reinhard-luminance |  |
| `L_{white}` | smallest scaled luminance mapped to pure white, eq. 4 (threshold in scaled luminance, same scale as L_sc) | 1 | `L_{white}` (colour.notation-reinhard-luminance) | colour.notation-reinhard-luminance |  |
| `x_{pix}` | pixel x coordinate in the image | 1 | `x` (colour.notation-reinhard-parameters); text/other: `(x,y)` | colour.notation-reinhard-parameters |  |
| `y_{pix}` | pixel y coordinate in the image | 1 | `y` (colour.notation-reinhard-parameters); text/other: `(x,y)` | colour.notation-reinhard-parameters |  |
| `N_{pix}` | total number of pixels in the image (paper symbol N) | 1 | `N_{pix}` (colour.notation-reinhard-parameters) | colour.notation-reinhard-parameters |  |
| `\delta_{L}` | small offset added to luminance to avoid the singularity at black pixels (paper symbol delta) | same as L_w | `\delta_{L}` (colour.notation-reinhard-parameters) | colour.notation-reinhard-parameters |  |
| `a_{key}` | key value, key value; target scaled luminance in eq. 2, before tone compression; paper symbol a; default 0.18 | 1 | `a_{key}` (colour.notation-reinhard-parameters) | colour.notation-reinhard-parameters |  |
| `X` | CIE 1931 tristimulus value X | same as Y | `X` (colour.notation-cie-xyz, colour.notation-stokes-cie-xyz) | colour.notation-cie-xyz |  |
| `Y` | CIE 1931 tristimulus value Y (luminance-like; photometric normalisation not stated in this source) | source-unspecified | `Y` (colour.notation-cie-xyz, colour.notation-stokes-cie-xyz) | colour.notation-cie-xyz |  |
| `Z` | CIE 1931 tristimulus value Z | same as Y | `Z` (colour.notation-cie-xyz, colour.notation-stokes-cie-xyz) | colour.notation-cie-xyz | Bare Z is never the ion charge (Z_ion). |
| `R_{sRGB}` | linear sRGB red tristimulus value | 1 | `R_{sRGB}` (colour.notation-srgb-linear) | colour.notation-srgb-linear |  |
| `G_{sRGB}` | linear sRGB green tristimulus value | 1 | `G_{sRGB}` (colour.notation-srgb-linear) | colour.notation-srgb-linear |  |
| `B_{sRGB}` | linear sRGB blue tristimulus value | 1 | `B_{sRGB}` (colour.notation-srgb-linear) | colour.notation-srgb-linear |  |
| `R'_{sRGB}` | non-linear (encoded) sRGB red value | 1 | `R'_{sRGB}` (colour.notation-srgb-nonlinear) | colour.notation-srgb-nonlinear |  |
| `G'_{sRGB}` | non-linear (encoded) sRGB green value | 1 | `G'_{sRGB}` (colour.notation-srgb-nonlinear) | colour.notation-srgb-nonlinear |  |
| `B'_{sRGB}` | non-linear (encoded) sRGB blue value | 1 | `B'_{sRGB}` (colour.notation-srgb-nonlinear) | colour.notation-srgb-nonlinear |  |
| `R_{8bit}` | 8-bit digital code value of the red channel | 1 (counts) | `R_{8bit}` (colour.notation-srgb-digital) | colour.notation-srgb-digital |  |
| `G_{8bit}` | 8-bit digital code value of the green channel | 1 (counts) | `G_{8bit}` (colour.notation-srgb-digital) | colour.notation-srgb-digital |  |
| `B_{8bit}` | 8-bit digital code value of the blue channel | 1 (counts) | `B_{8bit}` (colour.notation-srgb-digital) | colour.notation-srgb-digital |  |
| `\mathrm{WDC}` | white digital count | 1 (counts) | `WDC` (colour.notation-srgb-digital) | colour.notation-srgb-digital |  |
| `\mathrm{KDC}` | black digital count | 1 (counts) | `KDC` (colour.notation-srgb-digital) | colour.notation-srgb-digital |  |
| `x` | CIE chromaticity coordinate x | 1 | `x` (colour.notation-chromaticity) | colour.notation-chromaticity |  |
| `y` | CIE chromaticity coordinate y | 1 | `y` (colour.notation-chromaticity) | colour.notation-chromaticity | Never column depth (m) or pixel row (y_pix). |
| `z` | CIE chromaticity coordinate z | 1 | `z` (colour.notation-chromaticity) | colour.notation-chromaticity | Never redshift (z_g). |
| `\lambda` | wavelength column of the CMF file | nm (CMF table), cm in code | `wavelength` (colour.notation-cie-cmf) | colour.notation-cie-cmf | Wavelength only; affine parameter is lambda_aff, Compton length lambda_e. |
| `\bar{x}` | CIE 1931 2-deg colour matching function x-bar | source-unspecified (tabulated values) | `xbar` (colour.notation-cie-cmf) | colour.notation-cie-cmf |  |
| `\bar{y}` | CIE 1931 2-deg colour matching function y-bar | source-unspecified (tabulated values) | `ybar` (colour.notation-cie-cmf) | colour.notation-cie-cmf |  |
| `\bar{z}` | CIE 1931 2-deg colour matching function z-bar | source-unspecified (tabulated values) | `zbar` (colour.notation-cie-cmf) | colour.notation-cie-cmf |  |
| `K_{cd}` | luminous efficacy of monochromatic radiation of frequency 540e12 Hz, exactly 683 lm/W (SI defining constant) | lm W^-1 | `K_cd` (colour.k-cd-luminous-efficacy; BIPM SI Brochure, printed `Kcd`) | colour.k-cd-luminous-efficacy | (M1 2026-10-23) SI as printed (1 W = 1e7 erg s^-1). Subscript is a text label, braced per the spelling rules. |
| `Y_{pre}` | unexposed Y of the rendered spectrum (the Y of pipeline.display-colour before any exposure factor is applied) | same as Y | `Y_{pre}` (colour.display-exposure-convention) | colour.display-exposure-convention | (M1 2026-10-23) Replaces the bare L (and the retired synonym L_in) of the first display2 draft. Not the CIE tristimulus `Y` after exposure, not `Y_e`. |
| `t_{exp}` | exposure time factor of the stated camera model (logged per frame); t_exp > 0 | s | `t_{exp}` (colour.display-exposure-convention) | colour.display-exposure-convention | (M1 2026-10-23) Not `t`, `t_{ep}`, `t_{pr}`. |
| `k_{exp}` | exposure scale factor such that L_sc = t_exp k_exp Y_pre (logged per frame); k_exp > 0 | [Y_pre]^-1 s^-1 | `k_{exp}` (colour.display-exposure-convention) | colour.display-exposure-convention | (M1 2026-10-23) L_sc is dimensionless. The Reinhard key exposure is the special case k_exp t_exp = a_key/Lbar_w with t_exp = 1. |
| `I_{psf}` | specific intensity of the ideal Airy PSF at off-axis angle theta_psf | same as I | `I_{psf}` (rt.airy-intensity-pattern); `I(theta)` (source) | rt.airy-intensity-pattern | (M1 2026-10-23) |
| `I_{peak}` | central (theta_psf = 0) intensity of the Airy PSF | same as I | `I_{peak}` (rt.airy-intensity-pattern); `I_0` (source) | rt.airy-intensity-pattern | (M1 2026-10-23) `I_0` is the existing Bel02 surface-intensity row. Flux normalisation not harvested (entry notes). |
| `\theta_{psf}` | small off-axis angle in the camera frame, theta_psf = r_fp/f_foc (paraxial); distinct from the ray angle theta | rad | `\theta_{psf}` (rt.airy-intensity-pattern, rt.airy-first-zero, rt.psf-aperture-convention); `theta` (source) | rt.airy-intensity-pattern | (M1 2026-10-23) In the renderer, the offset from the centre of each source's PSF. Not `\theta_{field}` (pixel direction from the optical axis) or `\theta_{cam}`. |
| `x_{psf}` | Airy argument pi D_ap theta_psf/lambda | 1 | `x_{psf}` (rt.airy-intensity-pattern, rt.airy-first-zero, rt.psf-aperture-convention); `x`, `r'` (source, r' = pi D r/(f lambda)) | rt.airy-intensity-pattern | (M1 2026-10-23) Bare x is CIE chromaticity. |
| `r_{fp}` | radial distance in the focal plane | cm | `r_{fp}` (rt.airy-intensity-pattern); `r` (source) | rt.airy-intensity-pattern | (M1 2026-10-23) Not the radial coordinate `r`; `r_{pix}` is the same radius in pixels (r_pix = r_fp / pixel pitch). |
| `\mathrm{J}_1` | Bessel function of the first kind, order 1 | 1 | `\mathrm{J}_1` (rt.airy-intensity-pattern, rt.airy-first-zero, rt.psf-aperture-convention); `J_1` (source) | rt.airy-intensity-pattern | (M1 2026-10-23) Upright like `\mathrm{Ry}`: a special function, not a variable. Distinct from the Pot12 emissivity `J_{l}` (whose source spelling includes J_1) and `J`. |
| `D_{ap}` | circular aperture diameter of the camera model (logged parameter) | cm | `D_{ap}` (rt.airy-intensity-pattern, rt.airy-first-zero, rt.psf-aperture-convention); `D` (source) | rt.airy-intensity-pattern | (M1 2026-10-23) Bare `D` is distance. |
| `f_{foc}` | camera focal length | cm | `f_{foc}` (rt.airy-intensity-pattern); `f` (source) | rt.airy-intensity-pattern | (M1 2026-10-23) Physical focal length. The pixel-unit projection constant is `f_{cam}`. |
| `\theta_{1,psf}` | angular position of the first Airy zero, (x_1/pi) lambda/D_ap ~ 1.22 lambda/D_ap | rad | `\theta_{1,psf}` (rt.airy-first-zero); `1.22 lambda/D` (source, prose) | rt.airy-first-zero | (M1 2026-10-23) |
| `x_{1,psf}` | first positive root of J_1, the value of x_psf at the first Airy zero (3.8317059702 in rt.airy-first-zero) | 1 | `x_1` (rt.airy-first-zero) | rt.airy-first-zero | (M1 2026-10-23) Renamed from the proposed `x_1` to carry the psf scope of `x_{psf}` and `\theta_{1,psf}`; x_1 alone reads as a value of the CIE x (entry latex edit pending). |
| `K_{psf}` | unit-integral Airy kernel at wavelength lambda, K = A_psf / integral of A_psf over the stated support | sr^-1 | `K_\lambda` (rt.psf-aperture-convention) | rt.psf-aperture-convention | (M1 2026-10-23) Renamed from the proposed `K_\lambda`: subscript lambda means per unit wavelength here (`F_\lambda`, `B_\lambda`, `I_\lambda`); wavelength is an argument, K_psf(theta_psf; lambda). Entry latex edit pending. |
| `A_{psf}` | peak-normalised Airy pattern [2 J_1(x_psf)/x_psf]^2 at wavelength lambda (1 at the centre) | 1 | `A_\lambda` (rt.psf-aperture-convention) | rt.psf-aperture-convention | (M1 2026-10-23) Renamed from the proposed `A_\lambda`: per-unit-wavelength reading as for `K_{psf}`, and A_lambda is the usual spelling of extinction at wavelength lambda (cf. `A_V`). Entry latex edit pending. |

## Pipeline-level spellings not yet in notation entries

The `pipeline.*` entries and PRD prose were outside this pass, which edited notation entries only. They still use the source spellings below. Use the project symbol when you next touch them.

| pipeline spelling (entry) | project symbol | note |
|---|---|---|
| `y`, `T(y)` (atmosphere-solve, atmosphere-structure) | `m`, `T(m)` | column density |
| `u = r_g/R`, `1 + z`, `R_inf` (ns-exterior-spacetime, surface-parameters) | `x_g`, `1 + z_g`, `R^\infty` | |
| `R_inf`, `Ry`, `alpha` (constant-set) | `R_{Ry}`, `\mathrm{Ry}`, `\alpha_f` | Rydberg constant, Rydberg energy, fine-structure constant |
| `alpha_f`, `a-hat`, `q`, `m`, `b` (vacuum-polarisation) | `\alpha_f`, `\hat a`, `q_{vac}`, `m_{vac}`, `b_{rel}` | |
| `Theta` (surface-radiance, surface-radiance-h-atmosphere) | `\theta_{Bn}` | |
| `kappa^abs_j`, `kappa^sc_jj'` (mode-opacities) | `\kappa^a_j`, `\kappa^s_{jj'}` | |
| `beta (or q)`, `K_z,j`, `e_j,alpha` (polarisation-modes) | `\beta_{pol}`, `K_{z,j}`, `e_{j,\alpha}` | |
| `R_MM'`, `J_j = 1 - R_j` (surface-radiance-condensed) | `R_{jj'}`, `J_l = 1 - R_l` | mode vs linear-basis index |
| `g = nu_cam/nu_em` (frequency-shift) | `g_{shift} = \nu_{cam}/\nu` | bare g is surface gravity |
| `tau`, `E`, `L`, `u^mu`, `r_c`, `e_(a)` (rendered-frame, camera-worldline, camera-tetrad-projection, null-geodesic-trace) | `\tau_{cam}`, `\tilde E_{cam}`, `\tilde L_{cam}`, `u^\mu_{cam}`, `r_{cam}`, `e_{(a)}` | the camera's proper time, specific energy, angular momentum, four-velocity, radius and tetrad; tau, E and L are taken |
| `I_nu(x, y, nu)`, `Q, U, V`, `I_lambda`, `F_nu` | `I_\nu(x_{pix},y_{pix},\nu_{cam})`, `Q_\nu, U_\nu, V_\nu`, `I_\lambda`, `F_\nu` | spectral cube, Stokes parameters, flux |
| `f_col`, `p` (surface-radiance-blackbody) | `f_{col}`, `p_{ld}` | p belongs to the Pot12 fit exponents |
| `eta_i, eta_f`, `g-bar_ff(u, gamma^2)` (gaunt-factor) | `\eta_i, \eta_f`, `\bar g_{ff}(u_{kT},\gamma^2)` | Sommerfeld parameters; the subscript distinguishes them from the emissivity `\eta_\nu^j` |
| `B_pole`, `a` (key), `d` | `B_p`, `a_{key}`, `D` | |

## Remaining doubts

1. `\alpha` (Bel02 emission angle) and `\theta` (RT ray angle) are separate symbols even though they are equal for an emergent ray. Pech83 `delta` was merged into `\alpha`, and Pot12 `theta_k` into `\theta`.
2. `\nu` was given the meaning "local plasma frame". This absorbed the Bel02 `nu_0`, the PC03/McPHAC frequency and the frame-generic Pech83 `nu`. That is exact for a static star only. Rotation (Phase 4) will need separate static and comoving frequencies.
3. `\rho` merges the TOV gravitational mass density (energy/c^2) with the plasma rest-mass density. The two agree in the atmosphere and envelope but not in the core.
4. `\varepsilon_\perp`, `\varepsilon_\wedge` and `\varepsilon_\parallel` merge the vAL05 condensed-matter tensor with the PLCH04 plasma tensor. The sign of the gyrotropic element has not been checked across the two sources.
5. `k^\mu` merges two four-vector normalisations (Bel02 u_t = 1; Pech83 eq. 2.2). `\lambda_{aff}` likewise.
6. `T_{eff}^\infty` now includes Ho07's fitted T^infty for an emission-region atmosphere. `R_{em}^\infty` is kept apart from the whole-star `R^\infty`. `P_{spin}` merges the generic GH83 period with the RX J1856 value.
7. After r1, `F_{rad}` is identified with `F_{ph}` and `T_{eff,atm}` with `T_s` (VAL06 and McPHAC spellings kept as `sym_source`). Per-patch quantities that are only numerically equal are no longer kept as separate symbols.
8. `B_{dip}` (spin-down) and `B_p` (polar): the factor-2 convention between them has not been verified against a source.
9. `\lambda_e` is the reduced Compton wavelength hbar/(m_e c). The name follows PY01.
10. Pot12's emissivities `J`, `J_X`, `J_O` and `J_l` keep the letter J; `J_\nu` (with the nu subscript) is the mean intensity. They are visually close but no longer collide literally.
11. After r1 the compound sRGB rows, the PLCH04 copies of `E_V`, `\rho_V`, `E_1` and the doubled `\rho_1`, `\zeta_e` declarations were removed. Duplicate notation entries remain because no entry deletions were allowed. These are the two CIE-XYZ entries (colour.toml vs colour_srgb.toml), the compound and the split sRGB rows, and the vacuum symbols in both opacity.notation-vacuum and rt.notation-vacuum. They now carry identical `sym` and meaning.
12. Rejected notation entries were normalised too, because check.py includes them in its symbol table. Their adoption in this file is notational only and does not bless their content.
13. eos.notation-structure prints `Lambda` with PPP15's missing inverse. The canonical meaning states e^(-2 Lambda) = 1 - 2GM_r/(c^2 r), and the printed text is kept in `meaning_source`.
14. Non-notation entries (statements, latex, the pipeline prose) still use their source spellings, by design. The table above lists the pipeline cases; the equation entries keep the paper notation, with `sym_source` as the bridge.

## Review r1 outcome (truth/reviews/notation-r1.md)

The review returned a verdict for every dictionary row. This section records what was done and what is still open. `fix` rows whose remedy was a symbol, meaning, unit or alias edit were applied to the notation entries and this file. `fix` rows that need a different source, anchor or locator are marked "r1 fix pending" in the note column: those fields belong to the entry owners and were not edited. `reject` rows were kept and carry a "(disputed: ...)" note.

- **Rejected (kept, disputed):** `e` (reviewer's 2.7e-10 mismatch was at first not reproduced; re-verified 2026-10-23 with the exact hbar = h/(2 pi) it does reproduce and the reject is upheld as a limit on e^2 = alpha hbar c, see Conventions; the conversion row is kept as the project's definition), `\rho`, `\lambda_{aff}`, `k^\mu`, `\nu`, `T_{eff}^\infty`.
- **Applied:** source spellings and aliases (`h`, `A_{ion}`, `Z_{ion}`, `\varepsilon_{\perp,\wedge,\parallel}`, Pot12, GH83 and PPP15 App. B symbols); units (`rad` for Pech83 and GH83 angles, `W_L`, `L_{max}`, `L_{white}`, `u^\mu_{obs}`); meanings (`\mathrm{eV}`, `\mu_B`, `\Lambda`, `x_H`, `\rho_B`, `B_l`, `f_B`, `J_l` orientation, spectral Jacobians); compound rows split or removed; synonymous symbols merged (`\theta'_{sky}` into `\xi_{sky}`, `\Delta\phi_{orb}` into `\psi`, `T_{eff,atm}` into `T_s`, `F_{rad}` into `F_{ph}`, `I_X, I_O` into `I_\nu^j`); `K_{ij}` renamed `K_{ab}^{cond}`.
- **Pending (source or locator edits, not in this lane):** `N_A`, `c_1`, `c_2`, `\alpha_f`, `\gamma_B`, `f_{ion}`, `u_{kT}`, `\omega_{ci}`, `\zeta_p`, `\zeta_i`, `j`, `e_{j,\alpha}`, `q_{vac}`, `m_{vac}`, `\hat\kappa_\alpha`, `\nu^s_p`, `\nu_{p,\alpha}`, `\tilde\nu_\alpha`, `\gamma^2`, `T_s^\infty`, `T_b`, `g_{14}`, `B_{14}`.
- **Pending (structural):** merge of `\alpha` into `\theta`; fold of `I_\nu(m,\mu)`, `I_\nu(m,-\mu)` into `I_\nu`; BSk interior-fit rows (off the Milestones 1-3 render path) may be dropped from the active dictionary.

### Pipeline symbols with no canonical row (r1 cross-row finding)

Camera proper time, worldline coordinates and four-velocity, tetrad, specific energy and angular momentum with units, periastron and camera radius, photon four-momentum, direction vectors, frequency ratio, camera and emitter frequencies, Stokes `Q_\nu, U_\nu, V_\nu`, `I_\lambda`, `F_\nu`, surface field components, Sommerfeld parameters and unaveraged `g_{ff}`, parsec and solar-mass units (need a cited source), exposure, PSF and glare kernels. The earlier list under "Pipeline-level spellings" covers only part of this. No source anchors are available for these in the current notation entries, so no rows were invented.

### Open scientific doubts raised by r1

- The Pot12 linear-polarisation basis: the explicit vectors (App. B) put l = 1 perpendicular and l = 2 parallel to the plane of incidence, contradicting the paper's prose (checked against the local text).
- `\cos\psi_{vis}` and `b_{max}` hold only in the single-image regime; the validity limit R > 2 r_g was taken from the review and Bel02 and has not been re-derived here.
- Pechenick PDF pages 2 and 9 and PPP15 page 52 were image-checked by the reviewer only; scoped verification should be recorded before promoting their poor-OCR entries. The sRGB source labels itself obsolete.

## Additions 2026-10-23 (Milestone-1 harvest)

Merged from `truth/harvest/notation_additions_<lane>.md` for the lanes geodesic, camera, const, radiance, field_sky, sky, sky2 and display2, read as they stood at 13:10, after the camera, sky and display2 lanes revised their files following review. 102 rows were appended to the Dictionary, and each of their notes starts "(M1 2026-10-23)". By block: Constants 3, GR exterior 16, Camera frame and image plane (new) 24, Radiative transfer 1, Surface model 4, Astrometry and Gaia photometry (new) 38, Colour and display 16. Only `obs.notation-gaia-sky`, `obs.notation-gaia-photometry` and `obs.notation-gaia-propagation` are notation entries. For every other new row, the "adopted" column names the equation, constant or pipeline entry that introduces the symbol. The entry edits these decisions require are listed in `truth/harvest/notation_merge_report.md`. Until they are applied, `A_{u}` (obs.notation-gaia-sky) matches no row. The `\langle F_\lambda\rangle`, `\langle F_\nu\rangle` rename was already applied to obs.notation-gaia-photometry by the apply-sky-r2 lane, following the sky r2 review.

### Decisions

- `\tilde E`, `\tilde L` are adopted, not `E_{geo}`, `L_{geo}`. The subscript `geo` already means geometric units (`u_{geo}`). The reserved `\tilde E_{cam}`, `\tilde L_{cam}` are these symbols evaluated on the camera worldline.
- `\tau_{prop}` is the generic proper time and `\tau_{cam}` (reserved) is its value on the camera. Bare tau stays with the optical depths.
- `r_{peri}` is adopted, not the PRD's `r_p`, because subscript p already means proton and polar.
- `r_{cam}`, `e_{(a)}`, `u^\mu_{cam}`, `\nu_{cam}` and `g_{shift}` are the reserved Pipeline-level spellings, and now have rows.
- `\Theta_K(\theta)` is renamed `\Theta_K(\theta_{Sch})`. Its argument is the Boyer-Lindquist polar coordinate, and bare `\theta` is the RT ray angle.
- The Kerr rows (`a_{spin}`, `\Delta_K`, `\rho_K^2`, `\Sigma_K^2`, `\Theta_K(\theta_{Sch})`, `\mathcal{Q}`) are accepted only so that the a = 0 reduction can be read against DA09. The project sets a = 0 everywhere.
- The camera lane's static-observer label `st` is kept, although `u^\mu_{obs}` (Pech83, blessed) names the same observer. Under `obs`, `\nu_{obs}` would read as the observed (camera) frequency.
- `\alpha_{RA}` and `\delta_{Dec}` have one row each. The sky and field_sky lanes chose the same spelling but gave different meanings: Gaia ICRS at J2016.0 versus Walter+2010 with no stated epoch. The row meaning is generic, and the frame and epoch are attributes of each datum.
- `A_{u}` (sky) is renamed `\mathrm{au}` (const), so the astronomical unit has one spelling, written as a unit like `\mathrm{eV}`, `\mathrm{pc}` and `\mathrm{Ry}`. `A_{vel}` and `\tau_{A}` (sky2) are derived constants and keep their own rows.
- `\langle f_\lambda\rangle` and `\langle f_\nu\rangle` are renamed `\langle F_\lambda\rangle` and `\langle F_\nu\rangle`, as the sky-r1 reviewer also suggested. Flux density is `F_\lambda` (Kaplan row) and the reserved `F_\nu`, and lower-case `f_\nu` is the Eddington factor. The sky r2 review made the same decision, and it is applied in the entry.
- The display2 lane's revised spellings are accepted: `Y_{pre}`, `t_{exp}`, `k_{exp}`, and the psf-scoped `I_{psf}`, `I_{peak}`, `\theta_{psf}`, `x_{psf}`, `r_{fp}`, `\theta_{1,psf}`, `D_{ap}`, `f_{foc}`, `\mathrm{J}_1`. `L_{sc}` is reused: the display convention produces the Reinhard scaled luminance directly, and the key form is its special case.
- The rejected draft spellings are not added: `L_{in}` and `L_{abs}` from display2's first draft were superseded by `L_{sc}` and `Y_{pre}`. `I_0` would have clashed with the Bel02 row, and `\theta_1` would have read as a value of the RT ray angle.
- `x_1` is renamed `x_{1,psf}`, so that it carries the same psf scope as `x_{psf}` and `\theta_{1,psf}`.
- `K_\lambda` and `A_\lambda` are renamed `K_{psf}` and `A_{psf}`. In this dictionary a lambda subscript means per unit wavelength (`F_\lambda`, `B_\lambda`, `I_\lambda`), and A_lambda is the usual spelling of extinction (cf. `A_V`). Wavelength becomes an argument.
- `\mathrm{J}_1` is the Bessel function, because `J_1` is already a source spelling of the Pot12 emissivity `J_{l}`.
- `f_{cam}` (pixels, camera lane) and `f_{foc}` (cm, display2) are both kept. They are the same lens in different units: f_cam = f_foc / pixel pitch, and likewise r_pix = r_fp / pixel pitch.
- `R^{N}_\odot`, `L^{N}_\odot`, `q_{esu}` and `q_{C}` are rejected. Their entries (const.r-sun-nominal, const.l-sun-nominal, const.esu-from-coulomb) have status rejected. The first two were dropped because no Milestone 1-3 step uses them. The third was rejected because NIST prints only a rounded statC factor and does not identify the exact 10 c conversion with Gaussian charge.
- `\alpha_B` and `i` (field_sky) are rejected as rows. They are informal names for `\Psi_{pole}` and `\theta_{los}`, and bare `i` is a source spelling of the mode index `j`.
- The radiance call form `R^\infty(123\,\mathrm{pc})` is rejected. Ho07's radius is the emission-region radius `R_{em}^\infty` (collision row `R_inf`), as the lane's own later row says, so write `R_{em}^\infty(123\,\mathrm{pc})`.

### Collisions resolved (source spellings, this pass)

| source spelling | parties | resolution |
|---|---|---|
| `E` | rt.notation-val photon energy; DA09 specific energy (gr.timelike-*) | Photon energy keeps E. DA09 E -> `\tilde E`. |
| `L`, `L_z`, `L_in` | retired bare L; DA09 angular momentum; display2 first-draft luminances; IAU nominal L_sun | DA09 L_z -> `\tilde L`. Display -> `L_{sc}`, `Y_{pre}`. `L^{N}_\odot` rejected. |
| `R` | stellar radius; DA09 radial function R(r) | Stellar radius keeps R. DA09 -> `\mathcal{R}(r)`. |
| `tau` | optical depths; geodesic proper time; BL14 light time | `\tau_\nu`, `\tau_T` keep. -> `\tau_{prop}` (`\tau_{cam}` on the camera), `\tau_{lt}`, `\tau_{A}`. |
| `r`, `r_p`, `r_c`, `u_*` | radial coordinate (keeps r); Airy focal-plane radius; PRD periastron; pipeline camera radius; DA09 photon turning point | -> `r_{fp}`, `r_{peri}`, `r_{cam}`, `r_{tp}`. |
| `phi` | Bel02 surface azimuth (keeps `\phi`); DA09 coordinate; James phi_cs; surface-point azimuth | -> `\phi_{Sch}`, `\phi_{cam}`, `\phi_{surf}`. |
| `theta` | RT ray angle (keeps `\theta`); DA09 coordinate; James theta_cs; MG09 co-tetrad theta^(i); PY01 field angle; Airy off-axis angle; BL14 Delta theta | -> `\theta_{Sch}`, `\theta_{cam}`, `e^{(a)}`, `\theta_{Bn}`, `\theta_{psf}`, `\Delta\theta_{lt}`. |
| `alpha` | emission angle (keeps `\alpha`); Gaia and Walter+2010 RA; field_sky alpha_B; James lapse | -> `\alpha_{RA}`; alpha_B is `\Psi_{pole}`. The James lapse stays `alpha` only inside quoted James equations, with no row. |
| `delta` | Pech83 alias of the emission angle; Gaia and Walter+2010 Dec; DA09 delta_1 | -> `\delta_{Dec}`; `\delta_1` keeps its subscript. |
| `beta`, `gamma` | bending angle (keeps `\beta`); MG09/James camera speed; MG09 Lorentz factor | -> `\beta_{cam}`, `\gamma_{cam}` (`\gamma_B`, `\gamma^2` unchanged). |
| `b` | impact parameter (keeps b); James axial angular momentum p_phi (p_t = -1); BL14 position | James b: no row, distinct from b except for equatorial rays. BL14 -> `\vec{b}_{bary}`, `b_{bary}`. |
| `psi` | Bel02 orbital angle (keeps `\psi`); MG09 spatial-direction amplitude | MG09 psi: latex only, no row. |
| `n`, `B`, `N`, `n_F` | retired bare n; MG09 n and James B (direction of motion); James N; James n_F | -> `\hat v_{cam}`, `\hat N_{cam}`, `\hat n_{st}`. |
| `Lambda` | metric function (keeps `\Lambda`); MG09 boost matrix | -> `\Lambda_{boost}`. |
| `g` | surface gravity (keeps g); DA09 frequency ratio | -> `g_{shift}`. |
| `nu_0` | Bel02 emitter-frame frequency; DA09 observed frequency | Bel02 -> `\nu` (existing); DA09 -> `\nu_{cam}`. |
| `p^mu` | physical photon four-momentum (carries h nu); James normalised tangent (p_t = -1) | No row yet. They are related by p^(0) = E_obs/c = -p_mu u^mu/c (camera lane). |
| `f`, `f_c`, `<f_lambda>`, `<f_nu>` | Eddington factor `f_\nu`; James camera frequency; Airy focal length; Gaia mean fluxes | `f_\nu` keeps. -> `\nu_{cam}`, `f_{foc}` (pixel version `f_{cam}`), `\langle F_\lambda\rangle`, `\langle F_\nu\rangle`. |
| `D`, `d` | distance (keeps D); Airy aperture; James Doppler factor; Ho07 distance d | -> `D_{ap}`, `\mathcal{D}_{cam}`; Ho07 d is `D`. |
| `x`, `x_1` | CIE chromaticity (keeps x); PY01 compactness; WS83 x = r/r_g; Airy argument and its first zero | PY01 -> `x_g`; WS83 x stays in latex only (inverse of x_g at the surface); Airy -> `x_{psf}`, `x_{1,psf}`. |
| `a` | `a_{key}`, `a^{B.2}` etc.; DA09 Kerr parameter; PY01 dipole parameter | -> `a_{spin}`, `a_{dip}`. |
| `mu` | ray cosine (keeps `\mu`); WS83 dipole moment; BL14 total proper motion and its vector | -> `\mu_{mag}`, `\mu_{tot}`, `\vec{\mu}_{pm}`. |
| `chi` | PLCH04 polarisability (keeps `\chi`); PY01 magnetic colatitude; MG09 polar direction angle | PY01 -> `\theta_B`; MG09 latex only, no row. |
| `u` | retired bare u; BL14 source direction; MG09 four-velocity | -> `\vec{u}_{bary}`, `u^\mu_{cam}`. |
| `p`, `q` | Pot12 exponent and node index (keep); BL14 triad | -> `\vec{p}_{tri}`, `\vec{q}_{tri}` (and `\vec{r}_{tri}`). |
| `S`, `P`, `G`, `n_p` | retired bare S; pressure; Newton's G; proton density (keep); Gaia passband, pupil, G band, band flux | -> `S_{Gaia}`, `P_{pupil}`, `G_{Gaia}`, `n_{p,Gaia}`. |
| `A_V`, `A_u`, `A_lambda` | visual extinction (keeps `A_V`); BL14 au in km yr s^-1; Gaia au; Airy pattern | -> `A_{vel}`, `\mathrm{au}`, `A_{psf}` (kernel `K_{psf}`). |
| `varpi` | albedos `\varpi_\nu`, `\varpi_\nu^j` (keep); parallax | -> `\varpi_{plx}`. |
| `I_0`, `I(theta)` | Bel02 surface intensity (keeps `I_0`); Airy central and off-axis intensity | -> `I_{peak}`, `I_{psf}`. |
| `J_1` | Pot12 J_{l=1} (source spelling of `J_{l}`); Bessel function | -> `\mathrm{J}_1`. |
| `Delta`, `rho^2`, `Sigma^2`, `Theta`, `Q` | DA09 Kerr functions and Carter constant | -> `\Delta_K`, `\rho_K^2`, `\Sigma_K^2`, `\Theta_K(\theta_{Sch})`, `\mathcal{Q}`. |
| `xi` | retired bare xi; MG09 azimuth of a tetrad vector (gr.tetrad-vector-components latex) | Latex only, no row. |

### Alias bridges and notes on existing rows (no new row)

- `k_B`: Marr & Wilkin `k` (rt.planck-function-bnu, rt.planck-function-blambda).
- `\sigma_{SB}`: Marr & Wilkin `\sigma` (rt.planck-integral-stefan-boltzmann). The radiance lane's two sigma rows reduce to this one bridge.
- `D`: Ho07 `d` (surface.rxj1856-rinf-ho2007-rescaled-123pc). Ho07's scaled distance d_140 = D/(140 pc) has no row.
- `R_{em}^\infty`: Ho07 `R∞` (surface.rxj1856-rinf-ho2007-rescaled-123pc), with the call form `R_{em}^\infty(123\,\mathrm{pc})`.
- `\theta_B`: PY01 `\chi` (surface.py01-gr-dipole-magnitude, surface.field-colatitude-from-geometry).
- `\theta_{Bn}`: PY01 `\theta` (surface.py01-gr-dipole-field-angle). Identifying it with the field-normal angle is our step, which the field_sky lane checked against WS83 eqs. 12a-b.
- `x_g`: PY01 `x` (surface.py01-gr-dipole-a-of-x). WS83's x = r/(2M) = r/r_g is the inverse, equal to 1/x_g at the surface, so the two must not be mixed.
- `\lambda_{aff}`: DA09 `lambda` (eqs. 4, 8-10). For delta_1 = 1 it equals the proper time `\tau_{prop}` (gr.timelike-constants-schwarzschild). The r1 dispute stands.
- `GM/c^2`: DA09 `M (G=c=M=1)`.
- `r_g`: the exterior g_rr = (1 - r_g/r)^-1 is now carried by gr.schwarzschild-line-element (proposed).
- `\Lambda`: gr.schwarzschild-line-element (proposed) derives the exterior identity e^(-2 Lambda) = 1 - r_g/r (r >= R) from DA09 eq. 1 at a = 0. The PPP15 misprint stays in `meaning_source`.
- `\cos\psi`: Bel02 eq. 11 as used in gr.observed-flux-element-abs-jacobian, where psi is the total swept angle of an image ray and may exceed pi.
- `L_{sc}(x_{pix},y_{pix})`: also produced by colour.display-exposure-convention as L_sc = t_exp k_exp Y_pre. The Reinhard key form is the special case k_exp t_exp = a_key/Lbar_w.
- `\lambda`, `L_d`: display2 reuses them unchanged.

### Doubts opened by this pass

1. `st` and `obs` label the same local static observer (`e^{st}_{(t)}` = `u^\mu_{obs}`/c). Unifying them needs an edit to the blessed gr.notation-pechenick.
2. No local source prints the definition dtau^2 = -ds^2/c^2 behind `\tau_{prop}` (geodesic lane).
3. `K_{cd}` keeps SI units (lm W^-1), as the CODATA rows do. The Gaia rows are in CGS, with the printed SI units in the entries' `units_source`.
4. The blessed pipeline entry surface.field-colatitude-from-geometry writes `\chi`, `i` and `\alpha_B` in its own latex, which is not a paper spelling. The dictionary spellings are `\theta_B`, `\theta_{los}` and `\Psi_{pole}`.
5. Of the r1 cross-row list, these still have no row: the camera image-plane triad u_x, u_y (gr.camera-projection-convention, statement only), the photon four-momentum p^mu, the Stokes parameters `Q_\nu, U_\nu, V_\nu`, `I_\lambda`, `F_\nu`, the WS83 surface field components B^r-hat and B^theta-hat, the Sommerfeld parameters, a solar-mass unit, the pixel pitch, and the glare kernel.
