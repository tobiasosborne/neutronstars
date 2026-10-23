## harvest:colour skipped
- CIE 1931 colour-matching functions, CIE chromaticities of illuminants beyond D65 and 5778 K solar-white test: no CIE 015 source in refs/ (only refs/data/cie1931_2deg_cvrl.csv, which is not a section of this domain's assigned sources)
- IEC 61966-2-1 constants 3.2406255, 0.0031308, 0.04045: not in refs/; the v1 code used these without a local source
- Stokes eq. 1.3-1.6 (digital count scaling, 8-bit encoding): not needed for the HDR spectral-radiance pipeline; only 1.6 is in the equations transcription
- Reinhard eq. 5-9 (local dodging-and-burning operator, Gaussian kernels, phi sharpening): out of the assigned scope (global operator only)
- Reinhard physical units of L_w: the source states none, so the notation units are 'source-unspecified (luminance)'
- Stokes viewing-environment Table 0.1 (80 cd/m2, 64 lux, flare 1%, D50 ambient): out of scope; relevant to the PRD display model if reviewers want it

## harvest:eos skipped
- BSk24 coefficients: Potekhin et al. 2013 covers only BSk19, BSk20 and BSk21. BSk24 is not in any assigned source, so nothing was entered.
- Potekhin+2013 'Table 1' of M and R values: no such table exists (Table 1 lists Skyrme-force nuclear properties). The M-R benchmarks are text values in Sect. 6 and Fig. 8, entered as eos.bsk-*-benchmark.
- Ideal electron-ion plasma basics from PC01: refs/potekhin_chabrier_2001_eos.pdf is actually Potekhin & Yakovlev 2001 (A&A 374, 213). It has no degeneracy-parameter or Fermi-integral definitions beyond eq. 12 and Pr/Tr, which were entered. The real Potekhin & Chabrier 2001 EoS paper is not in refs/.
- Tables 3 and 4 (b_i, c_i) for BSk19 and BSk20. BSk21 b_i entered; BSk21 c_i and Table 5 (d_i) not entered.
- Eq. 8 (rho(H) fit) with Table 5, and Eq. 9 (Ye, Ymu core fractions) with Table 6: not needed for Milestones 1-3.
- Eqs. 10-21 and Tables 7-10 (crust composition, nuclear shapes, conductivity parameters): out of scope.
- Table 1 nuclear-matter properties (n_caus etc.): not needed.
- PPP15 thermal-structure content (eq. 37 Tb-Ts, eq. 38, Appendix B eqs. B.1-B.5 Tp/Teq fits) and PY01 eqs. 10-11 (GR dipole field), eqs. 27-29 and the Greenstein-Hartke Ts(theta) relation: belong to the surface domain, not eos.
- Separate dM_r/dr = 4 pi r^2 rho entry: PPP15 prints only the integral definition of M_r, so it is recorded inside eos.enclosed-mass-and-rg and eos.tov-pressure.
- Kinematic notation for the plasma Coulomb parameter Gamma: the fit's adiabatic index Gamma is noted as unrelated, but no plasma-Gamma entry exists in these sources.

## harvest:gr skipped
- SOURCE PROBLEM: refs/riley_2023_xpsi.pdf/.txt is the wrong paper. It is arXiv:2308.14312, 'Local-Global Pseudo-label Correction for Source-free Domain Adaptive Medical Image Segmentation', not X-PSI (the JOSS paper), so nothing was harvested from it. refs/INDEX.md and its checksum need fixing and the correct X-PSI paper needs fetching.
- Timelike geodesic constants of motion (energy and angular momentum for free fall): not stated in any of the 3 usable sources.
- Static-observer orthonormal tetrad and proper time: only the timelike leg u^mu (Pechenick eq. 2.4) is stated. The spatial legs and dtau = e^Phi dt are not written explicitly, so they are not entered.
- z and 1+z: no source defines them. Beloborodov and Pechenick only give nu/nu_0 = (1 - r_g/R)^(1/2). The 1+z convention is noted in gr.frequency-redshift and is not a notation symbol.
- Spectral intensity transformation I_nu(nu) = (1 - r_g/R)^(3/2) I_nu0(nu_0): it follows from combining gr.invariant-inu-over-nu3 with gr.frequency-redshift, but no source states it. It should go in a derivation layer, not as a sourced entry.
- Exterior lapse e^{2Phi} = 1 - r_g/r for r > R: Potekhin gives only the surface boundary condition, and Pechenick gives no explicit line element. gr.metric-static-spherical notes this.
- Beloborodov eq. 2 (series for y = 1 - cos psi, x^2 term vanishes), eq. 3 (approximate trajectory r(psi)), eq. 4 (cosine-approximate flux), and Sect. 4 pulsar formulae (eqs. 5-8): not needed for Milestones 1-3, or derivable from entries already present. Eq. 4 follows from eq. 11 plus eq. 1.
- Pechenick eqs. 2.1 (angle between two photons), 2.3 (radial photon four-velocity), 2.7 (orbit equation in r; equivalent to 2.11 and covered by gr.photon-sphere), 3.x polar-cap light-curve machinery, 4.1, and Appendix A2-A15 and A18-A26 (tidal cone distortion, flash-sphere timing): outside scope or redundant.
- Potekhin eq. 39 (T_eff defined as an area integral) and the TOV equations (eq. 3): they belong to the surface and eos domains; surface.ppp15-teff-definition and eos.tov-pressure already exist.
- Dexter and Agol (Kerr) and Doppler boosts: Phase 4, and not in the four assigned sources.

## harvest:opacity skipped
- Explicit cold-plasma dielectric elements epsilon, g, eta (PLCH04 only cites Ginzburg 1970; the explicit forms are van_adelsberg_lai_2006 eqs. 10-12, not in the assigned sources and built on Meszaros 1992 Gaunt factors). The tensor structure is recorded but the elements are not, so the opacity domain cannot yet compute mode vectors from a plasma state.
- PC03 Appendix B Born formulae B20-B22 for Lambda^e_alpha and the quantum correction factor for eq. 51 (out of scope: 'not every appendix detail'). Eq. 49 (normalised Lambda^e) skipped for the same reason.
- PC03 eqs. 46-48 (proton-proton absorption and the Lambda_pp fit) and eq. 54: PC03 says pp is negligible, so only an approximation entry is recorded. v1 implemented eqs. 47-48 anyway.
- PC03 eq. 45 (customary ion free-free recipe) is called erroneous by the paper; mentioned in notes only.
- PC03 eqs. 25-26 (Gnedin-Pavlov q, p parameters and the |q|>>|p| approximation): superseded by the PLCH04 K_j, beta parametrisation. Eq. 26 appears only in a note.
- PC03 eqs. 29-30 (diffusion-approximation effective opacities): belongs to the RT domain.
- PC03 eqs. 35-36 and bound-bound/bound-free cross sections: partial ionisation is a Phase 4 item.
- PLCH04 eqs. 7-17 (oscillator-strength sums, Kramers-Kronig, perturbative atomic model): needed only for partially ionised plasma.
- Vacuum resonance adiabatic-conversion and P_jump formulae (vAL06 eqs. 3-4, SPW09 eqs. 16-17): mode-conversion belongs to RT.
- McPHAC eqs. 16-17 (Thomson dipole and isotropic phase functions): scattering redistribution is RT-domain. Mentioned in notes only.
- Mészáros (1992) ground-Landau-level opacity formulae, Compton correction and the (nu/nu_c)^4 X-mode scaling: Meszaros 1992 is not in refs/, so nothing from it was recorded.
- Lai (2001) hbar omega_cp = 6.3 B12 eV, SPW09 eq. 12 E_ci = 0.635 (B/1e14 G)(Z/A) keV and PHC16 zeta_i = 0.0737 (Z/A) B12/T6: dropped to stay within 40 entries; the proton-cyclotron relations omega_cp = (m_e/m_p) omega_ce and beta_p ~ 0.0732 B12/T6 are in notation notes.
- Numerical values of sigma_T and fundamental constants: belong to the const domain.

## harvest:rt skipped
- Planck function B_nu and its Rayleigh-Jeans limit: none of the three sources states the formula (they only say B_nu is the Planck/blackbody intensity). It belongs in const/colour from a source that prints it.
- Stefan-Boltzmann value and F = sigma T_eff^4 as a standalone constant: no numeric value in these sources (SPW09 eq. 18 and VAL06 eq. 9 use it and are entered as flux-constraint equations).
- Limb-darkening law: no source gives a definition or a fitted law. Only the emergent I(0,mu) = 2P(0,mu) is entered (rt.emergent-intensity-mcphac).
- Grey / Hopf starting solution: SPW09 says only 'a grey temperature distribution' (prose); Hopf function is not in any of the three sources. Only the Eddington profile with T(0) = 0.265 T_eff is entered.
- Eddington-factor isotropic values f = 1/3, h = 1/2 (McPHAC Sect. 4.5): own entry would have broken the 40-entry cap. They are mentioned in the notes of rt.eddington-factors-mcphac but are not covered by that entry's anchor, so a reviewer should treat them as unanchored until promoted.
- McPHAC eq. 11 (reduced Thomson opacity), eq. 12 (free-free opacity) and the Gaunt factor: belong to the opacity lane.
- McPHAC unmagnetised validity bound B << 1e10 (T/1e6 K) G, no-self-irradiation (R > 1.5 R_S) and no-electron-conduction assumptions: not entered, cap reached.
- McPHAC discretisation details A7-A11 (finite differences) and boundary rows A26, A28-A33: only A12-A23 core is entered. The T_eff-dependent grid choices (K=100, M=40, tau_max > 80, 1e-4 convergence) are numerical settings, not physics.
- SPW09 eq. 11 (ideal-gas EOS), eqs. 12-13 (ion cyclotron energy, equal-energy grid), and the opacity formulae in VAL06 (eqs. 10-29: dielectric tensor, mode polarisation vectors, ff and scattering opacities): opacity lane.
- VAL06 Stokes transfer eq. 39-41 (matrix M and emission vector): matrix text is garbled in the text layer; only the mode-Stokes conversion relations (eqs. 43-47) are entered. Polarisation propagation and the observed-polarisation sections (eqs. 49-87) are Milestone-5 or out of PRD scope.
- VAL06 eq. 37-38 (finite-difference form of the transfer equation for Runge-Kutta): numerical scheme, not physics.
- Kurucz (1970) details (Lambda_diag(tau) dependence, Avrett-Krook and surface corrections) and Ho & Lai 2003 f_B (eq. 2.41): not in refs/; cited by SPW09 and VAL06 but cannot be harvested.

## harvest:surface skipped
- Dipole field components B_r, B_theta and |B| vs magnetic colatitude: none of the three assigned sources states them. GH83 eq. 1 is the polar-cap radius, not a field magnitude. PHC16 only cites a 'relativistic dipole [72]' (Ginzburg & Ozernoi) with no formula. Needs another source in refs/.
- Magnetic colatitude of a general surface point from obliquity and rotation phase: only GH83 eq. 21 (pole-to-line-of-sight angle) is stated, and it is included. The surface-point dot-product formula in archive/v1 dipole.jl is not in any assigned source.
- Cosine-squared ansatz T_eq + (T_pole - T_eq) cos^2(theta_B): not in GH83 (see notes). The v1 claim that GH83 give T ~ |cos theta|^(1/4) is also not in the paper. The combined closed form of eqs. 14+15 (K_perp=0 gives T^4 = T_pole^4 * 4cos^2(beta)/(1+3cos^2(beta))) is my own algebra, not a source statement, so I left it out. A reviewer can derive it as a consistency check.
- Pulsed fraction as a general definition: no assigned source defines it. GH83 eq. 23 degree of modulation M = (peak-min)/min is included as surface.gh83-modulation-definition, and it is not the (max-min)/(max+min) convention. Suleimanov 2009 and Ho 2007 mention pulsed fractions but were not assigned.
- GH83 eqs. 16-20 and 24 (photon count rate, J(N_H,T), Table 1): instrument/ISM-specific Einstein-IPC machinery, not needed for Milestones 1-3.
- GH83 eq. 25-26 (Radhakrishnan-Cooke position angle), PSR 1055-52 and PSR 1509-58 fits: radio-pulsar geometry and pulsar-specific fits, out of scope.
- PPP15 eq. 38 and the rest of the envelope structure equations, neutrino emissivities, conductivities: heat-transport microphysics, not the surface map.
- PHC16 light-bending approximation cos(theta_k) = x_g + (1-x_g) cos(theta) (cites Beloborodov [70]): belongs to the gr domain; left to the GR harvester.
- PHC16 spindown field B_s ~ 1.5e13 G for RX J1856 and R_inf = 17.2 km, z_g ~ 0.25: obs domain, not harvested here.
- Pulsed-fraction upper limits for RX J1856 (Ho 2007, Burwitz 2003): sources not assigned to this domain.

## review-r1:colour
- Silent collision: Reinhard's scaled luminance `L(x,y)` conflicts with stellar luminosity `L` in eos.toml:29. Pixel `(x,y)` also shares symbols with chromaticity coordinates without a scope declaration.
- Logarithm conventions need explicit scoping: BSk defines `log=log10` at refs/potekhin_2013_bsk_eos.txt:267-268; do not silently apply that convention to Reinhard's `exp/log` expression.
- The equation sidecar silently repairs both power-law images at refs/stokes_1996_srgb_equations.txt:17-23. It is not faithful corroborating evidence.
- Accepted matrices do not establish the absolute spectral-radiance-to-XYZ conversion or its normalisation. That missing bridge remains necessary for the first scientifically grounded colour image.

Orchestrator notes. (a) The reviewer's claim that srgb13.gif and srgb21.gif repeat the red input in every row is false. I read both images directly, and the G and B rows use G and B (G' and B' for 1.7b). So colour.srgb-transfer-power stays proposed with the verdict "fix" recorded as "disputed". (b) colour.srgb-inverse-transfer-power is rejected. Its statement wrongly claims a G/B copy-paste typo, and it silently repairs the printed unprimed branch condition (srgb20.gif prints R_sRGB,G_sRGB,B_sRGB > 0.03928). A replacement entry is needed that transcribes the printed condition. (c) colour.srgb-inverse-transfer-linear was fixed to record the printed unprimed condition, with the primed reading flagged as editorial. It is still proposed. (d) All 8 blessed entries with ocr=poor got image_verified = true. The reviewer reported vision passes on the page 4 image and the GIFs, and on the HTML tables for the two tables. (e) The reinhard-log-average-luminance 1/N-outside-exp printed form is blessed as a flagged transcription, so the reboot must not silently pick the geometric-mean form.

## review-r1:eos
INCOMPLETE. The reviewer for eos.toml had not finished, so truth/reviews/eos-r1.md does not exist yet. Only the progress log truth/reviews/eos-r1.md.log is there. No verdicts were applied and truth/entries/eos.toml is unchanged, so all 33 entries are still status=proposed. I did not run `python3 -I truth/check.py` this step, so check_errors=-1 means not run. To finish, wait for the codex process to exit and confirm that truth/reviews/eos-r1.md exists. Then apply the per-entry verdicts and run the check.

## review-r1:gr
Reviewer cross-entry issues (verbatim, plus my notes). (1) Mechanical check: initial execution exited 0 and printed `318 entries: {'proposed': 318}; 0 errors, 215 warnings`. Warnings comprised 212 missing provenance classes and the declared `alpha`, `lambda`, `phi` conflicts. (2) Shared schema blockers: all 35 non-notation GR entries lack `provenance`; all 23 poor-OCR entries lack `image_verified`. Relevant PDF pages were visually inspected: Beloborodov 4; Pechenick 1-4, 9, 11; Potekhin 4, 41. Scientific acceptance above does not resolve these metadata omissions. (3) Metric consistency: using the printed metric and Pechenick's radial null vector (line 122), at 1-r_g/R=2/3, expected g(v,v)=0; computed -0.833333333, error 0.833333333, tolerance 1e-12. The erroneous radial factor (Potekhin e^{2Lambda}=1-2GM_r/(c^2 r), which should be (1-2GM_r/(c^2 r))^-1) also appears in `eos.notation-structure`. (4) Numerical source discrepancy: two substitutions of Pechenick equation 2.13 agree on the limb values -0.887573189,-0.635279431,-0.485939428. Absolute errors against Beloborodov are 0.001573189,0.002279431,0.001939428, all exceeding 0.0005. His quoted visible fraction 0.742 similarly differs from computed 0.742969714 by 0.000969714. (5) Hidden notation collisions: glyph/string differences conceal emission-angle `alpha` versus fine-structure alpha and photometric alpha; photon u^k versus observer u^mu; potential B versus magnetic field B; photon versus spin nu; affine parameter versus wavelength lambda. Beloborodov's covariant constants also need an explicit convention bridge to Potekhin's metric. (6) Silent synonyms: reconcile `Phi`/`Phi(r)`, `Lambda`/`Lambda(r)`, and the radius aliases `R_inf`/`R^inf`/`R^infty`; distinguish apparent radius from the Rydberg constant R_infty. (7) Need: no clear drops. These claims support Milestones 1-3 as physics or verification targets. Published numerical benchmarks require declared recomputation routes and tolerances; infinity-observer formulas cannot implement the finite moving camera (PRD.md:30). My notes: (a) I independently integrated Beloborodov eq. 9 and confirmed mu_v = -0.88757, -0.63528, -0.48594 for R/r_g = 2, 2.5, 3 (visible fraction 0.74297), so the two rejects (gr.metric-static-spherical, gr.visible-mu-exact-values) are justified; the reject of the metric entry needs a replacement entry with the correct g_rr = (1-2GM_r/(c^2 r))^-1 from a supportable source (Pechenick/Misner-style), and eos.notation-structure should be re-reviewed for the same Potekhin e^{2Lambda} typo. (b) The codex_review.sh log (gr-r1.md.log) contained a '-s: command not found' error from the shell (script was likely being edited concurrently by another agent), yet gr-r1.md was written at 07:27:37 with gr-specific content, so I used it; rerun if the provenance of that run matters. (c) Reviewer flagged 'ocr=poor' accepted entries as PDF-checked; I set image_verified = true on those accepted/fixed poor entries where the reviewer said 'PDF checked'. (d) Provenance fields remain missing on all non-notation gr entries (check.py only warns); blessed gr entries still carry the provenance and ownership warnings. (e) alpha/lambda/phi conflict warnings between gr.notation-geometry and gr.notation-pechenick remain for the notation consolidation pass; gr.notation-pechenick alpha was narrowed to the polar-cap meaning and notes the Appendix cone half-angle use.

## review-r1:opacity
- Notation collisions: opacity `theta_B` means ray-field angle; surface `θ_B` means magnetic colatitude. Bare `h` in the McPHAC formula also conflicts with EOS enthalpy; bare ν conflicts with observational spin frequency. Canonical mappings must survive ASCII/Greek spelling differences.
- Units and aliases: `e`, `m_e`, `m_p` use Gaussian-cgs units here versus SI in `const`. Register conversions and aliases such as `omega_ce`/`omega_c`, `alpha_f`/`α`, and `kappa_nu`/`κ_ν`; the checker misses these.
- Mass convention: PC03 defines m_H = m_p + m_e, whereas McPHAC calls m_H the hydrogen mass and keeps it separate in Eq. 11. Do not silently identify these conventions (PC03:38, McPHAC:227-232).
- Coordinate contract: Eq. 18 vectors use z ∥ k; opacity cyclic components use z ∥ B. The rotation remains an explicit implementation prerequisite (PLCH04:343-353, PC03:279-285,382-383).
- Dependency cycle: `mode-scattering-opacity` and `approx-isotropic-scattering` depend on each other.
- All reviewed subject matter can support Milestone 3; no clearly unneeded entry warrants `drop`.
- Five notation bundles (notation-fields, notation-particles-plasma, notation-photon-cyclotron, notation-vacuum, notation-opacity-symbols) were rejected because the single anchors do not support the bundles. They must be re-split into independently anchored definitions (B_14, gamma_B, b_rel, a_m; Z/A, Y_e, alpha_f, f_ion; McPHAC nu and PHC16 zeta_e, zeta_i; E_V, rho_V, f_B, E_1; kappa-hat_alpha, sigma_T,r, gamma^2) in the consolidation pass.
- Partial dispute: for approx-ff-proton-collisions-negligible the reviewer cited ratio nu^pp/nu^ff,p, but PC03 Eq. 54 reads nu^{ff,p}/nu_alpha and nu^pp/nu^s_p, so source wording was used and beta_p < 1 added.
- All 10 ocr=poor entries are among the fix verdicts (all stay proposed). The reviewer stated it read the page image for each, and image_verified = true was set on them, but none are blessed yet.
- Mechanical check: 0 errors for the whole truth/ tree. Warnings come from other domains: missing provenance/ownership, and gr alpha/lambda/phi symbol conflicts.

## review-r1:rt
- Canonical notation remains unresolved: RT μ is a positive half-range variable in H/V, signed in S, and differs from GR μ = cos ψ. Photon ν also collides with timing ν; McPHAC j labels angular points while magnetic-transfer j labels modes.
- Merge documented aliases explicitly: y/y_col/m; g/g_s/g_surf; θkB/αBk/opacity θB; EBi/Eci; and the duplicated vacuum symbols. Greek, ASCII and LaTeX spellings currently conceal matches.
- V:186 delegates the actual function fB to Ho & Lai (2003), absent locally. These entries establish neither a computable fB nor permission to substitute fB = 1.
- All 39 entries have a plausible Milestone 1–3 role; no clear `drop`.
(Executor note: rt.feautrier-matrix-coeffs-mcphac and rt.source-function-val were rejected; the McPHAC matrix form and VAL06 general scattering source function therefore have no blessed entry. Fixed-but-proposed entries need a round 2 accept before blessing. No ocr=poor entry was accepted; te-rybicki got image_verified=true but remains proposed.)

## review-r1:surface
- Appendix-B polar `B_12` conflicts semantically with local-field `B_12` in `eos` and `opacity`. Adopt a distinct project symbol such as `B_p,12`, preserving source notation in LaTeX. (Applied in surface.notation-envelope-fit: sym renamed to "B_p,12"; the schema check had flagged the B_12 conflict across eos.notation-electron-gas, opacity.notation-fields and the surface entry.)
- Consolidate apparent-radius aliases `R^∞`, `R^inf`, `R_inf`, and `R_\infty`; explicitly map `σ_SB`, `σ_sb`, and `σ`, including their SI/CGS units.
- Poor-OCR equations were checked against PDF images in this review; record `image_verified = true` before blessing (done for the 11 accepted ocr=poor entries). The current schema also requires provenance classifications for all 34 non-notation entries, with reconstruction routes for retained fits (currently only warnings: 'provenance' and 'ownership' fields missing).
- Reviewer's first mechanical check reported a 'Y' symbol conflict between colour.notation-cie-xyz and colour.notation-stokes-cie-xyz (concurrent changes by other files; the later recheck passed).
- surface.notation-envelope-fit: g_14 = g/(10^14 cm s^-2) is defined explicitly only in refs/potekhin_chabrier_2001_eos.txt L540, not in PPP15; cross-source symbol.

## review-r2:colour
Reviewer verbatim: (1) "Pixel `(x,y)` still needs explicit scoping against chromaticity coordinates (truth/entries/colour_srgb.toml:101). XYZ/RGB declarations also duplicate those in colour_srgb.toml:10." (2) "All five reviewed non-notation entries require `provenance` and `ownership` before blessing. Both image-based entries lack `image_verified=true`; their local GIFs received a vision pass here." (3) Mechanical check noted three GR notation conflicts: alpha, lambda, phi (gr.notation-geometry vs gr.notation-pechenick). Actions taken: added provenance=convention/ownership=leaf to the accepted non-notation entries; set image_verified=true on colour.srgb-transfer-power (reviewer vision-read srgb12/13.gif). Rejected: colour.srgb-inverse-transfer-linear (its statement says 0.0393 to 4 s.f.; 12.92*0.00304 = 0.0392768, which is 0.03928) - needs a corrected statement and a new review. Not done: the colour.notation-cie-xyz vs colour_srgb.toml duplication and the (x,y) scoping are left for the notation consolidation pass.

## review-r2:gr
Reviewer, verbatim:
- Rejected (surviving content): gr.notation-pechenick: Appendix PDF p. 854 shows Greek alpha, not Latin a; OCR misleading at Pechenick:520-526. Surviving: alpha denotes polar-cap angular radius in Sections II-III and is reused for photon-cone half-angle in the Appendix.
- gr.observed-flux-element: "Per image branch" overclaims the signed derivative expression: secondary branches can reverse its sign. Surviving: printed equality on the primary monotonic branch; general branches require |d cos(alpha)/d mu|. Counterexample: integrating Pechenick eq. 2.12 at R/M=3.2, alpha=85deg gives psi=190.698757308deg. With I_0 dS/D^2=1, expected absolute-Jacobian flux +0.031152015859; signed expression -0.031152015859; error 0.062304031717, tolerance 1e-10.
- gr.approx-observer-at-infinity: Locator p. 846 omits eq. 2.6 on p. 847. The proposed finite-radius replacement is derived for a static observer, not the moving camera. Surviving: stationary-infinity approximation; eq. 2.5 gives the finite-radius angle only in the static frame.
- Notation: unresolved collisions include GR alpha versus atomic alpha, and photon u^k versus observer u^mu. Silent synonyms remain Phi/Phi(r), Lambda/Lambda(r), and R_inf/R^inf/R^infinity.
- Schema: all 14 reviewed non-notation entries lack provenance and ownership (the 14 accepted or rejected non-notation entries; check.py currently reports 0 errors, so this is only a warning-level gap). Record this review's image verification for gr.notation-pechenick and gr.approx-opaque-surface; both lacked image_verified (I set image_verified = true on gr.approx-opaque-surface since it is blessed; gr.notation-pechenick is rejected).
- Dependencies: gr.metric-static-spherical remains rejected; dependent accepts do not rehabilitate it or resolve the component-convention bridge. No clear drops.
- check.py warnings still include declared conflicts for alpha, lambda, phi between gr.notation-geometry and gr.notation-pechenick (524 warnings total, 0 errors).
Note: all 3 rejected entries (gr.notation-pechenick, gr.observed-flux-element, gr.approx-observer-at-infinity) remain in gr.toml with status rejected and a round-2 reject review; no round 3 per instructions. Poor-OCR accepted entries (null-geodesic-constants, null-geodesic-components, max-deflection-integral, disk-angular-radius, invariant-inu-over-nu3, approx-schwarzschild-exterior, approx-opaque-surface) all had image_verified = true and reviewer said PDF confirmed (schwarzschild-exterior note cites text-line matches only, but image_verified was already true).

## review-r2:opacity
- Initial mechanical run, exit 0: `357 entries: {'blessed': 137, 'rejected': 28, 'proposed': 192}; 0 errors, 532 warnings`. Warnings include missing provenance/ownership and conflicting GR symbols `alpha`, `lambda`, `phi`.
- Register the synonym opacity `theta_B` ≡ RT `θ_kB` ≡ `α_Bk`; distinguish it from surface `θ_B` (magnetic colatitude) and `θ_Bn` (field–normal angle). The checker's literal symbol strings conceal these collisions.
- `mode-scattering-opacity` and `approx-isotropic-scattering` depend on each other, violating the required dependency DAG.
- Several accepted equations still depend on rejected notation bundles; their source transcription does not repair those dependencies.

Rejected in round 2 (no round 3), with surviving content: notation-modes-polarization (cyclic components need frame designation z'||B; K_j uses z||k), approx-ground-landau-level (drop 'any field strength' overreach), ff-effective-frequency (reduction to McPHAC Eq.12 needs nu_e,0 << omega), approx-ff-proton-collisions-negligible (Eq.54 ratio is nu^pp/nu_alpha^ff,p), approx-isotropic-scattering (McPHAC flux claim needs separate anchored entry), dielectric-tensor-polarizability (add dipole approximation), vacuum-coeffs-fit (expansions A7-A9 present locally in Ho-Lai L214-229), vacuum-resonance-energy (f_B defined locally in Ho-Lai L444-457), ff-opacity-kramers-gaunt-mcphac (equality with PC03 Eqs 37+41 is a weak-damping limit). All ocr=poor accepted entries cite page-image confirmation in the reviewer's reason. Final check: 357 entries, blessed 211, rejected 44, proposed 102 (file-wide; opacity r2 subset all resolved), 0 errors.

## review-r2:rt
Incomplete: the codex round-2 review of rt.toml was started in the background but had not finished when this was forced to return. truth/reviews/rt-r2.md did not exist yet. No verdicts were read and truth/entries/rt.toml was not modified. All 15 targeted ids remain proposed. truth/check.py was not re-run, so check_errors is not a measured value. The rerun needs to be done.

## review-r2:surface
INCOMPLETE: the round-2 review had not finished when this output was forced. No output file exists yet, only surface-r2.md.log. No verdicts were read and truth/entries/surface.toml was not modified. All 15 round-2 ids are still proposed. check.py was not run after the review (-1 means not run). Re-run this step, or read truth/reviews/surface-r2.md once codex finishes.