# GAPS: work order for the next harvest sessions

Date: 2026-10-23. State: `python3 -I truth/check.py --render --all` gives 387 entries
(277 blessed, 80 rejected, 30 proposed), 0 errors and 3 warnings. The tree has 146 green,
230 yellow and 11 red nodes. The 30 proposed entries are the whole `pipeline.*` spine. It has
not been reviewed: `truth/reviews/pipeline-r1.md` does not exist, only its `.log`.

Inputs: PRD.md, DEPENDENCY_AUDIT.md (Parts B and C), harvest/spine_gaps.md,
harvest/wf1_cross_entry_issues.md, harvest/wf2_cross_entry_issues.md, wfA_results.json,
wfB_results.json, TREE.md, INDEX.md and refs/INDEX.md.

**Sources on disk.** There are 51 keys. The 17 recomputation-chain papers fetched today are:
`vanhoof_2014_gaunt`, `karzas_latter_1961_gaunt`, `hummer_1988_gaunt`,
`potekhin_1994_coulomb_landau`, `heyl_hernquist_1997_qed`, `kohri_yamada_2002_vacuum`,
`lai_ho_2002_vacuum_resonance`, `ho_lai_2003_atmospheres_vacuum`,
`potekhin_1999_magnetic_conductivity`, `baiko_yakovlev_1995_ephonon`,
`slattery_doolen_dewitt_1982_ocp`, `wilms_allen_mccray_2000_ism`, `verner_1996_photoionization`,
`morrison_mccammon_1983_ism`, `cardelli_1989_extinction`, `predehl_schmitt_1995_nh_av` and
`wasserman_shapiro_1983_gr_dipole`. None of the 17 has a truth entry yet.

**Ranking.** Gaps are ordered by milestone. Within a milestone they are ordered by how many
spine nodes the gap blocks: the node that needs it plus every spine node above it.
`pipeline.constant-set` blocks 29 nodes and `pipeline.rendered-frame` blocks 1.

**Column key.**
- "Local" is a refs key plus a section or equation, or a tex label (`l.` is a line in `.tex`
  or `.txt`).
- "Fetch" gives author, year, venue and ID. "(unv.)" marks an ID that was not checked against
  a file.
- "Est." is the estimated number of new or replacement entries.

---

## 1. Blocks Milestone 1 (spacetime, camera, blackbody disc, real sky)

| # | Missing | Spine node (blocks) | Local source / fetch | Est. |
|---|---------|---------------------|----------------------|------|
| 1.1 | **Constants:** IAU exact au and pc, and IAU 2015 nominal GM_sun (datum). **SI to Gaussian-CGS bridge**, e.g. e[esu] = e[C]·c/10 (derivation). Today the CODATA entries are SI only. | constant-set (29) | **Fetch:** IAU 2012 Res. B2 (au), IAU 2015 Res. B2 (pc) and B3 (nominal solar values), from iau.org PDFs; Prša et al. 2016, AJ 152, 41, arXiv:1605.09788 (unv.). The Gaussian-unit appendix is in Jackson, *Classical Electrodynamics* 3rd ed. (1999, ISBN 0-471-30932-X), which is a manual download. | 5 |
| 1.2 | **Exterior Schwarzschild line element with the correct g_rr = (1−r_g/r)^−1.** All metric entries are rejected because of the PPP15 e^{2Λ} misprint. A replacement is needed for the parents of `gr.frequency-redshift` and `gr.approx-schwarzschild-exterior`. | ns-exterior-spacetime (10) | **Local:** `dexter_agol_2009_kerr` §2, eq. 1 (Boyer–Lindquist Kerr line element; add an a=0 reduction as a derivation entry). Cross-check: `wasserman_shapiro_1983_gr_dipole` §II (Schwarzschild background, y=2M/r, l.195–209, OCR poor). | 2 |
| 1.3 | **Timelike geodesics:** conserved E and L, radial potential, periastron condition R(r_p)=0, proper-time quadrature. | camera-worldline (9) | **Local:** `dexter_agol_2009_kerr` §2, eqs. 3–6: R(r) and Θ(θ) with δ₁=1 for timelike orbits, the t and φ quadratures, and λ(r). Add an a=0 derivation entry. The periastron and E>1 hyperbolic conditions are derivable. | 4 |
| 1.4 | **Full static orthonormal tetrad** (only u^μ exists today), **local Lorentz boost and relativistic aberration**, and the **camera projection** (pinhole or equidistant). | camera-tetrad-projection (8) | **Local:** `wasserman_shapiro_1983_gr_dipole` l.258 (Schwarzschild orthonormal tetrad, OCR poor, needs an image check). **Fetch:** boost and aberration from James, von Tunzelmann, Franklin & Thorne 2015, CQG 32, 065001, arXiv:1502.03808 (unv.), or Rybicki & Lightman 1979, *Radiative Processes*, Ch. 4 (book, ISBN 0-471-82759-2). The projection is a convention: a pipeline entry with source PRD, which needs a PRD decision. | 5 |
| 1.5 | **Dipole field:** B_r, B_θ, \|B\|, the field–normal angle Θ (tan Θ = ½ tan θ_B in flat space), the **GR dipole**, and the magnetic colatitude of a surface point from obliquity and phase (A28, A30). | surface-maps-field (7) | **Local:** `potekhin_yakovlev_2001_cooling` eq. labelled `dipole` (tex l.479–491): B(χ)=B_p√(cos²χ+a²sin²χ), tan θ = a tan χ, a(x). The flat limit a→½ gives the v1 formula. Cross-check: `wasserman_shapiro_1983_gr_dipole` eqs. 12a–b. The colatitude follows from spherical trigonometry (derivation, source PRD). | 4 |
| 1.6 | **General frequency ratio** g = (p·u_cam)/(p·u_em), and the **special-relativistic Doppler factor** for the R2 far-field check. | frequency-shift (6) | **Fetch:** Rybicki & Lightman 1979 §4.1–4.9 (book), or James+2015 (as 1.4). g itself is derivable from blessed `gr.invariant-inu-over-nu3`. | 2 |
| 1.7 | **Bending integral from a finite camera radius r_c**, and the **multi-image flux Jacobian \|d cos α/dμ\|** to replace the rejected `gr.observed-flux-element`. | null-geodesic-trace (6) | **Local derivation:** from blessed `gr.null-geodesic-radial` and `pechenick_1983` eq. 2.12. The counterexample in wf1 review-r2:gr (R/M=3.2, α=85°, ψ=190.70°) becomes a verification datum. | 3 |
| 1.8 | **Planck function B_ν(T) as an entry**, with the analytic check ∫B_ν dν = σT⁴/π. | surface-radiance-blackbody (6) | **Not printed locally.** CODATA gives only c₁L and c₂; SPW09 and LH02 l.438 only name B_ν. **Fetch:** Rybicki & Lightman 1979, eq. 1.51 (book). | 2 |
| 1.9 | **Camera PSF / diffraction model** (Airy pattern for a stated aperture), used for point-source stars and the unresolved neutron star. | star-scale-adaptive (5), sky-radiance | **Fetch:** Born & Wolf, *Principles of Optics* §8.5.2 (book, manual). | 2 |
| 1.10 | **M1 parameter labelling.** Ho+2007 T∞, R∞ and z are code output at d=140 pc; add a derivation entry rescaling R∞ ∝ d to 123 pc. State the alternative input: the Chandra blackbody kT∞ and R∞/d (both blessed). | surface-parameters (3) | **Local:** `ho_2007_rxj1856_thin_atm` Table 1, l.39–41; `walter_2010_rxj1856_distance`; `burwitz_2003`. | 1 |
| 1.11 | **Sky.** RX J1856 RA and Dec. Gaia DR3 catalogue (positions, parallaxes, G/BP/RP, XP spectra) and its astrometric model. Passbands and zero points. Bright stars missing from Gaia (G ≲ 3). Starless diffuse Milky Way map and its licence. | sky-radiance (3) | **Local:** `walter_2010_rxj1856_distance` tex l.480 (α=18h56m35.11s, δ=−37°54′30.5″, "nominal", epoch unstated). **Fetch:** Gaia Collab., Vallenari+2023, A&A 674, A1, arXiv:2208.00211; Riello+2021, A&A 649, A3, arXiv:2012.01916 (passbands, zero points); Montegriffo+2023, A&A 674, A3, arXiv:2206.06205 (calibrated XP spectra, measured, which avoids a synthetic-SED table); Lindegren+2021, A&A 649, A2, arXiv:2012.03380 (astrometric solution); van Leeuwen 2007, A&A 474, 653, arXiv:0708.1752 (Hipparcos, bright stars); diffuse-map candidates Mellinger 2009, PASP 121, 1180, arXiv:0908.4360, and Leinert+1998, A&AS 127, 1 (all unv.). The catalogue rows themselves are measured data with checksums. | 10 |
| 1.12 | **Display colour.** A blessed **inverse sRGB power segment** (both versions are rejected; re-harvest with the printed unprimed condition). The tristimulus integral X = k∫S x̄ dλ and its normalisation. **K_cd = 683 lm/W** for absolute radiance to luminance. | display-colour (3) | **Local:** `stokes_1996_srgb` eq. 1.7a/b (srgb20.gif, srgb21.gif; condition printed with unprimed symbols, threshold 0.03928). **Fetch:** BIPM *SI Brochure* 9th ed. (2019), free, defines K_cd; tristimulus normalisation from CIE 015:2004 (paywalled, manual) or Wyszecki & Stiles 1982 (book). | 4 |
| 1.13 | **Exposure model** (physical exposure H ∝ L·t/N², or a key-based model) and the **inverse Reinhard operator** (derivable). Veiling glare and video time remapping are PRD M4 and are listed here only so they are not lost. | display-exposure-tone (2) | Inverse operator: derivation from blessed `colour.reinhard-simple-operator`. **Fetch:** exposure from ISO 12232 (paywalled) or a photographic-optics text (book); glare from Spencer, Shirley, Zimmerman & Greenberg 1995, SIGGRAPH (unv.), at M4. | 3 |

Polarisation-basis transport is not listed here. M1 radiance is an unpolarised blackbody with
Q=U=V=0, so the item is deferred to 3.5.

---

## 2. Milestone 2 (Earth-distance photometry, R5; own parameter fit)

| # | Missing | Spine node (blocks) | Local source / fetch | Est. |
|---|---------|---------------------|----------------------|------|
| 2.1 | **Measured per-band HST optical/UV magnitudes of RX J1856** with photometric systems. Locally there is only Walter F475W (system unconfirmed). Kaplan 2011's 27.13 is a prediction and Table 5 holds fitted summaries. Keep the Ho 0.12 mag and Kaplan 0.05 mag extinctions separate (obs r1). | earth-photometry (4) | **Fetch:** van Kerkwijk & Kulkarni 2001, A&A 380, 221, astro-ph/0110065 (unv.), the source Kaplan 2011 uses for RX J1856; Walter & Matthews 1997, Nature 389, 358 (unv.). | 6 |
| 2.2 | **HST filter throughput curves** for the bands in 2.1 (WFPC2, STIS, ACS F475W): measured data files. | earth-photometry (4) | **Fetch:** STScI synphot/CDBS throughput tables (no paper; checksum each file). STMAG is local (`obs.stmag-definition`). | 4 |
| 2.3 | **Interstellar extinction:** the CCM89 a(x), b(x) curve and the A_V(N_H) relation, both fits, harvested as stubs. Own route: a grain model plus Mie theory. | earth-photometry (4) | **Local:** `cardelli_1989_extinction` §3, eqs. 1–3 (OCR poor); `predehl_schmitt_1995_nh_av` N_H/A_V = 1.79e21 (l.51, l.1050; OCR poor). **Fetch (own route):** Weingartner & Draine 2001, ApJ 548, 296, astro-ph/0008146; Draine 2003, ARA&A 41, 241, astro-ph/0304489 (unv.); Bohren & Huffman 1983 (book, Mie). | 5 (+deep) |
| 2.4 | **Chandra LETG spectra and responses** (ARF/RMF): measured data. | surface-parameters (3) | **Local:** ObsIDs 113, 3380, 3381, 3382, 3399 in `burwitz_2003_rxj1856_xray` Table 1 (tex l.140–158). **Fetch:** the data from the Chandra archive, checksummed. | 3 |
| 2.5 | **X-ray ISM photoabsorption σ_ISM(E)** (A34 L2). | surface-parameters (3) | **Local:** `wilms_allen_mccray_2000_ism` §2.1–2.4 (gas, molecules and grains; abundances); `verner_1996_photoionization` §3, eq. 1 + Table 1 (fit, stub); `morrison_mccammon_1983_ism` Table 2 (verification only). Own route through Opacity Project photoionisation is deep and deferred: Seaton+1994, MNRAS 266, 805 (unv.). | 8 |
| 2.6 | **Likelihood and fit procedure** (Poisson C-statistic for grating counts, χ² for photometry), and the pulse-amplitude forward model. `obs.rxj1856-pulse-amplitude` is blessed. | surface-parameters (3) | **Fetch:** Cash 1979, ApJ 228, 939 (unv.). The pulse-amplitude model is a derivation over pipeline nodes (source PRD). | 3 |

---

## 3. Milestone 3 (surface models; pass R5 with one model)

| # | Missing | Spine node (blocks) | Local source / fetch | Est. |
|---|---------|---------------------|----------------------|------|
| 3.1 | **Gaunt factor (A4, first M3 blocker).** Exact g_ff(η_i,η_f) from Sommerfeld's hypergeometric form; Maxwell average ḡ_ff(u,γ²); the non-Born ratio ḡ/ḡ_Born. The gffgu.dat verification fixture needs an entry (table, verification only, columns swapped). Re-harvest the B=0 Kramers form. | gaunt-factor (13) | **Local:** `vanhoof_2014_gaunt` §2.1–2.2 (tex l.129–326, exact calculation), §3 (l.524, thermal average), App. A (series); `karzas_latter_1961_gaunt` §II, eqs. for free–free (OCR poor, needs an image check); `hummer_1988_gaunt` (fit, verification only); `refs/code/McPHAC/gffgu.dat` plus `haakonsen_2012` footnote 2 (column swap); re-harvest `haakonsen_2012` eq. 12 with the ν_e,0 ≪ ω qualifier. **Optional fetch:** Bethe & Salpeter 1957 (book), for an independent ODE check of the Coulomb functions. Note: the audit gives arXiv:1407.1083 for van Hoof; the file on disk is 1407.5048. | 12 |
| 3.2 | **Magnetic free–free corrections (A16, A18).** PC03 App. A (Λ_pp chain A1–A6 and eq. 48). App. B Coulomb matrix elements B6–B8 and finite-proton-mass Λ^e_α B19–B22. The §4.3.2 quantum correction factor. Landau-state matrix elements. Landau wavefunctions. | free-free-opacity (12) | **Local:** `potekhin_chabrier_2003_ff_opacity` App. A, App. B, eqs. 46–48 and 54, §4.3.2; `potekhin_1994_coulomb_landau` Appendix, eqs. A3–A10 (l.701+, OCR poor); Landau functions from `potekhin_1994` §2 or `ruder_1994_atoms_strong_B` (OCR degraded). Re-harvest `opacity.ff-effective-frequency` (weak-damping qualifier) and `opacity.approx-ff-proton-collisions-negligible` (ratio ν^pp/ν_α^{ff,p}). | 16 |
| 3.3 | **Magnetic scattering redistribution** (mode j→j′, angle-dependent; A9). Break the cycle `opacity.mode-scattering-opacity` ↔ `opacity.approx-isotropic-scattering`. For Compton/Kompaneets, record only an approximation entry (kT/m_ec² ~ 1e-4 for RX J1856). | scattering-opacity (12), mode-opacities (9) | **Fetch:** Ho & Lai 2001, MNRAS 327, 1081, astro-ph/0104199 (unv.; paper I of the local HL03); Mészáros 1992, *High-Energy Radiation from Magnetized Neutron Stars* (book, ISBN 0-226-52094-3, manual). Coherent scattering: local `rt.approx-coherent-scattering-mcphac` (blessed); Kompaneets from Rybicki & Lightman Ch. 7 (book). | 6 |
| 3.4 | **Vacuum polarisation (A20).** One-loop Heyl–Hernquist closed form as an independent check of the PLCH04 A2–A6 entries. Weak- and strong-field expansions. **A derivation of the constant 0.0329199** (PLCH04 tex l.1047). Kohri–Yamada as a numerical check. | vacuum-polarisation (12) | **Local:** `heyl_hernquist_1997_qed` §2 (l.99–455) and §2.1 series and asymptotics; `ho_lai_2003_atmospheres_vacuum` eqs. labelled `vpalo`/`vpqlo`/`vpmlo` and `vpahi`/`vpqhi`/`vpmhi` (l.460–480); `kohri_yamada_2002_vacuum` §2–3, fitting formula `fitting_mu2_2` (l.407), Fig. 2. Re-harvest `opacity.vacuum-coeffs-fit` without the false "absent" note. **Optional fetch:** Euler & Heisenberg 1936, Z. Phys. 98, 714; English translation arXiv:physics/0605038 (unv.). | 9 |
| 3.5 | **Plasma dielectric tensor and normal modes (A19)**, plus **polarisation transport** (moved here from M1). Explicit cold-plasma elements ε, g, η (the opacity harvest flagged them as missing). A tensor entry that states the dipole approximation. Frame-explicit z′∥B versus z∥k. Parallel transport of the polarisation basis along Schwarzschild rays. | polarisation-modes (11) | **Local:** `ho_lai_2003_atmospheres_vacuum` §2.1, eqs. `epsij0`, `epsii0`, `epszz0`, `epsxy0` (l.386–405) and the mode parameter `polarb0`/`polarbvp` (l.596–624); `potekhin_lai_chabrier_ho_2004_polarization` eqs. 1–5 (re-harvest `opacity.dielectric-tensor-polarizability` with the dipole approximation). **Fetch** (transport): Poutanen 2020, A&A 641, A166, arXiv:2006.08287 (unv.); Pineault 1977, MNRAS 179, 691 (unv.). | 8 |
| 3.6 | **Atmosphere EOS.** Ideal fully ionised H, ρ = m_p P/(2k_BT) (A7). An m_H versus m_p convention entry (PC03 versus McPHAC). Validity flags for Landau quantisation and non-ideality. The Fermi-integral and μ closure that PY01 eq. 12 omits. | atmosphere-structure (10) | **Local:** `suleimanov_potekhin_werner_2009_mag_atm` eq. 11 (ideal gas, tex l.369; skipped by the rt and opacity lanes); `potekhin_chabrier_2013_mag_eos` (magnetised, fully ionised); `potekhin_chabrier_2000_eos` (non-magnetic, with OCP fit). The audit's "PC01" points to a deleted mislabelled file; these two local papers replace it. **Fetch** if the closure is needed: Chabrier & Potekhin 1998, PRE 58, 4941 (unv.). | 5 |
| 3.7 | **Vacuum resonance and mode conversion (A21–A23).** f_B. E_ad and γ_res. A **Landau–Zener derivation**. The **two-mode amplitude evolution equations** (L2 check). | vacuum-mode-conversion (9) | **Local:** `ho_lai_2003_atmospheres_vacuum` §2.3 eq. 2.41 f_B (txt l.444–457), E_V eqs. 2.38–2.40, §2.4 eqs. `adiabat`/`enadiabat` and K± (eq. 2.43); `lai_ho_2002_vacuum_resonance` §3 (γ_res, P_jump, Δz_res, l.246–340); `van_adelsberg_lai_2006` §3 (conversion inside the RT). Re-harvest the rejected `opacity.vacuum-resonance-energy` and `rt.vacuum-resonance-val`; their only defect was a false note that HL03 is absent. **Fetch:** Zener 1932, Proc. R. Soc. A 137, 696, or Haxton 1995 (as cited by LH02) (unv.); Lai & Ho 2003, ApJ 588, 962, astro-ph/0211315 (mode evolution, unv.). | 8 |
| 3.8 | **Feautrier solver.** An implementable McPHAC coupling matrix (re-harvest eqs. 13–15 with the identity term and the conversion to the dμ′ kernel). A replacement for the defective VAL06 eq. 30 source function (use SPW09). Elimination and boundary discretisation. The exact Hopf function (grey limit). The grey Eddington seed as the A12 L0 step. Unsöld–Lucy correction. | atmosphere-solve (8) | **Local:** `haakonsen_2012_mcphac` eqs. 10 and 13–15; `suleimanov_potekhin_werner_2009_mag_atm` source function (tex l.313–356); `potekhin_2012_condensed_surface` eqs. 31–34 (Unsöld–Lucy, skipped by the condensed lane). **Fetch:** Mihalas 1978, *Stellar Atmospheres* (book, manual); Hopf function from Mihalas or Chandrasekhar 1960, *Radiative Transfer* (Dover, ISBN 0-486-60590-6). | 8 |
| 3.9 | **Envelope T_b→T_s chain (A29, about 1200 LoC)**, **or** a parametrised T_s map fitted to the pulse amplitude (Part C item 7; this needs a decision). Chain: PPP15 eqs. 26, 31–36, 38, A.1–A.3; magnetised Coulomb logarithms; electron–phonon scattering; OCP MC with S(q) and free energy; magnetised EOS; Fe/C bound-free opacity (avoid with an H/He envelope). | envelope-tb-ts (7), surface-maps-temperature (6) | **Local:** `potekhin_pons_page_2015_thermal` (tex labels `kappa-eff`, `omegatau`, `sigma-tau`, `tau-quant`, `nu_ei_nondeg`, `nu_ee_nondeg`, `tau`, eq. 38); `potekhin_1999_magnetic_conductivity` §3–4 and fits `taupar`/`tauperp`/`Lpar`/`Lperp` (tex l.1005–1105); `baiko_yakovlev_1995_ephonon` §4 fits (l.437); `slattery_doolen_dewitt_1982_ocp` Table I (verification); `potekhin_chabrier_2000_eos` (OCP fit, Debye–Hückel limit); `potekhin_yakovlev_2001_cooling` (envelope fits, verification). **Fetch:** an OCP MC/Ewald method from Frenkel & Smit (book), or Brush, Sahlin & Teller 1966, J. Chem. Phys. 45, 2102 (unv.); Opacity Project (Seaton+1994, unv.) only if the envelope is Fe/C. A parametrised map is one derivation entry with source PRD. | 25 (or 2) |
| 3.10 | **Condensed surface.** General Fresnel reflectivity at arbitrary θ_k and θ_B (vAL05 eq. 13, App. A/B; P12 App. A quartic A.6–A.21, root selection A.23 stated for n_z,j, boundary conditions A.27–A.29). Melting/condensation criterion. Damping rates from Potekhin 1999. Surface density from Medin & Lai. | surface-radiance-condensed (6) | **Local:** `van_adelsberg_2005_condensed_surface` §2.1–2.3, eq. 13, App. A/B; `potekhin_2012_condensed_surface` App. A–B; `potekhin_1999_magnetic_conductivity` (τ_∥, τ_⊥). Five rejected rt.cond entries must be re-harvested (see 4.1). **Fetch:** Medin & Lai 2006, PRA 74, 062507, astro-ph/0607166; Medin & Lai 2007, MNRAS 382, 1833, arXiv:0708.3863 (unv.). | 18 |
| 3.11 | **H-atmosphere emergent intensity and verification fixtures.** Re-harvest I(0,μ) = 2P(0,μ) with the corrected locator. McPHAC spectra as a fixture. SPW09 and HL03 figure values. | surface-radiance-h-atmosphere (6) | **Local:** `haakonsen_2012_mcphac` p. 5, §3.3 (the unnumbered equation before eq. 17); `refs/code/McPHAC` runs; `suleimanov_potekhin_werner_2009` figures 4–6 and `ho_lai_2003` §5 (digitised, `table` provenance, verification only). | 6 |
| 3.12 | **Decision: partial ionisation.** The spine assumes fully ionised H, which CLAUDE.md says is valid for T_eff > 5e5 K. RX J1856 has T∞ ≈ 4.3e5 K, and the Ho+2007 fit is partially ionised at B ~ 4e12 G. R5 may need bound–bound and bound–free opacities. | surface-radiance-h-atmosphere (6), condensed (6) | **Local:** `potekhin_chabrier_2003` (partially ionised EOS and opacities), PLCH04 eqs. 7–17, `potekhin_1994` (atom in B). **Fetch:** Potekhin & Pavlov 1997, ApJ 483, 414 (photoionisation, unv.). Needs a PRD decision before harvest. | 1 decision (+~20) |

---

## 4. Database hygiene

### 4.1 Rejected entries that must be re-harvested (still needed)

"Live deps" means non-rejected entries or spine nodes that still depend on the rejected entry.

| # | Rejected id | Defect, and what survives | Live deps | Action |
|---|-------------|---------------------------|-----------|--------|
| 1 | `gr.metric-static-spherical` | PPP15 prints g_rr wrongly | `gr.frequency-redshift`, `gr.approx-schwarzschild-exterior` | Replace via 1.2 and repoint the deps. |
| 2 | `eos.notation-structure` | Same g_rr typo; bundle not covered by its anchor | `eos.surface-gravity` | Split it. Source g_14 from PY01 l.540. |
| 3 | `opacity.ff-opacity-kramers-gaunt-mcphac` | Equality with PC03 needs ν_e,0 ≪ ω | spine free-free-opacity | Re-propose with the qualifier (3.1). |
| 4 | `opacity.ff-effective-frequency` | The McPHAC reduction is unqualified | `opacity.damping-frequencies`, `opacity.ff-cross-section-*`, spine | Re-propose (3.2). |
| 5 | `opacity.approx-ff-proton-collisions-negligible` | Eq. 54 ratio misquoted | spine free-free-opacity | Re-propose with ν^pp/ν_α^{ff,p} and β_p < 1. |
| 6 | `opacity.dielectric-tensor-polarizability` | Dipole approximation not stated | `opacity.mode-parameters-k`, `opacity.mode-parameter-beta`, `opacity.approx-vacuum-low-energy`, spine | Re-propose (3.5). |
| 7 | `opacity.vacuum-resonance-energy` | False note that HL03 is absent | `opacity.vacuum-resonance-density` | Re-propose and link HL03 eq. 2.41. |
| 8 | `rt.vacuum-resonance-val` | Same false note | `rt.adiabatic-energy-val`, `rt.vacuum-relevance-val` | Re-propose. |
| 9 | `opacity.vacuum-coeffs-fit` | False "expansions absent" note (they are in HL03 l.214–229 txt) | none (verification target) | Re-propose as a verification fit. |
| 10 | `opacity.approx-isotropic-scattering` | McPHAC flux claim not separately anchored; forms a cycle | `opacity.mode-scattering-opacity` | Split it and break the cycle. |
| 11 | `opacity.approx-ground-landau-level` | "Any field strength" overreach | none | Re-propose, narrowed. |
| 12–17 | `opacity.notation-fields`, `-particles-plasma`, `-photon-cyclotron`, `-modes-polarization`, `-vacuum`, `-opacity-symbols` | Bundles not covered by a single anchor | 14 blessed opacity entries | Re-split into about 15 independently anchored notation entries (list in wf1 review-r1:opacity). |
| 18 | `rt.emergent-intensity-mcphac` | Wrong locator | spine h-atmosphere | Re-propose (3.11). |
| 19 | `rt.feautrier-matrix-coeffs-mcphac` | Identity term and kernel unresolved | none | Re-propose (3.8). |
| 20 | `rt.source-function-val` | VAL06 eq. 30 mixed indices (a source defect) | none | Replace with the SPW09 form (3.8). Keep VAL06 as a "disputed" note. |
| 21 | `rt.cond-emissivity-definition` | Polarisation basis contradiction (P12 §2.2.2 versus App. B) | `rt.cond-reflectivity-from-amplitudes`, `-linear-polarisation-degree`, `-inner-bc` | Resolve the basis, then re-propose. |
| 22 | `rt.cond-damping-rates` | False "Potekhin 1999 absent" note | `rt.cond-approx-damping-low-energy` | Re-propose with `potekhin_1999` links. |
| 23 | `rt.cond-fit-coefficient-data` | Ambiguous `1/2 sqrt(1+B13)` | spine condensed | Re-propose with 1/(2√(1+B13)); the 1/3 multiplies sin²θ_B. |
| 24 | `rt.notation-cond-energies-geometry` | Locator stops before p. 5 | 4 rt.cond entries | Extend the locator to §2.3.1, p. 5. |
| 25 | `rt.notation-cond-emissivity-fit` | W_L units, basis labels, per-keV Jacobian | `rt.cond-kirchhoff-law`, `rt.cond-inner-bc` | Re-propose. |
| 26 | `surface.notation-envelope-fit` | Cites the deleted PC01 txt | 5 surface entries | Re-source g_14 to PY01 txt l.540. |
| 27 | `colour.notation-stokes-cie-xyz` | "Luminance-like" Y not supported | `colour.srgb-colorimetric-basis`, `colour.srgb-stokes-xyz-to-linear-matrix` | Merge with `colour.notation-cie-xyz` and fix both together. |
| 28 | `colour.notation-reinhard-luminance` | L_d ∈ [0,1) is false for eq. 4; stilb claim unsupported | `colour.reinhard-log-average-luminance`, `colour.reinhard-scaled-luminance` | Re-propose, scoped to eq. 3. |
| 29 | `colour.srgb-inverse-power` / `colour.srgb-inverse-transfer-power` (duplicates) | Printed unprimed condition not stated | spine display-colour (gap) | One replacement (1.12). Drop the duplicate. |
| 30 | `colour.srgb-inverse-transfer-linear` | Threshold 0.0393 should be 0.03928 | none | Re-propose. |
| 31 | `obs.approx-vkk2008-dipole-spindown` | Sign of Ė with ν̇ < 0 | `obs.rxj1856-spin-period` (wrongly), `obs.rxj1856-dipole-field` | Remove the spin-period dependency. Re-propose with an explicit sign convention. |
| 32 | `gr.observed-flux-element` | Signed Jacobian on secondary branches | none (spine avoids it) | Re-propose with the absolute value (1.7). |
| 33 | `gr.approx-observer-at-infinity` | Locator misses eq. 2.6; static frame only | none | Re-propose as a static-frame statement. |
| 34 | `gr.visible-mu-exact-values` | Printed values fail their rounding | none | Keep rejected. Add our computed −0.887573189, −0.635279431, −0.485939428 (visible fraction 0.742969714) as a derivation/verification datum. |
| 35 | `const.notation-fundamental`, `gr.notation-pechenick` | Anchor and symbol mismatch | none | Re-split, low priority. |

### 4.2 Disputed anchors and source conflicts

| # | Item | Finding | Action |
|---|------|---------|--------|
| 36 | `eos.bsk-mmax-central-density`: is the BSk20 Mmax 2.27 or 2.16 M_sun? | **No error in the entry.** The anchor spans two clauses: "ρc = 2.69 × 1015 g cm−3 for BSk20, and Mmax = 2.27 M⊙ at". The 2.27 belongs to BSk21. tex l.1509–1512 gives BSk20 Mmax = 2.16 at ρc = 2.69e15 and BSk21 Mmax = 2.27 at ρc = 2.27e15, consistent with the values [3.48, 2.69, 2.27]e15. The wfB tex lane misread the anchor. The entry is dropped (A24), so the issue is moot for M1–M3. | Close the dispute. If the entry is ever revived, use an anchor that names a single EOS and add anchor_tex `$\rhoc=2.69\times10^{15}$`. |
| 37 | D65 white y = 0.3291 (stokes l.621) versus 0.3290 (l.773) | Unresolved within the source | Record both. Pick one in a convention entry and state which. |
| 38 | `colour.srgb-transfer-power` | Blessed, but the reviewer verdict "fix" is recorded as "disputed" | A third-party image check of srgb12/13.gif. |
| 39 | Beloborodov 2002 limb values and visible fraction 0.742 | Off by 1.5–2.3e-3 from exact integration | Covered by 34. |
| 40 | P12 polarisation basis (§2.2.2 versus Fig. 1 and App. B) | Source contradiction | One decision entry. Blocks 21 and 25. |
| 41 | PC03 m_H = m_p + m_e versus McPHAC m_H | Two conventions | A notation entry for each, plus a mapping (3.6). |
| 42 | PC03 venue | CLAUDE.md says A&A 399, 1007; the local file is ApJ 585, 955 | Fix CLAUDE.md. |
| 43 | van Hoof arXiv id | The audit says 1407.1083; the file on disk is 1407.5048 | Fix DEPENDENCY_AUDIT.md. |
| 44 | Audit "PC01" (Potekhin & Chabrier 2001 EOS) | The file was a misfiled copy of PY01 and has been deleted | Re-point the audit to `potekhin_chabrier_2000_eos` and `potekhin_chabrier_2013_mag_eos`. |

### 4.3 Graph and schema defects

| # | Defect | Action |
|---|--------|--------|
| 45 | Cycle `opacity.mode-scattering-opacity` ↔ `opacity.approx-isotropic-scattering`. It is the cause of the 11 red nodes, including the root. | Fix through 4.1 row 10. |
| 46 | 25 rejected entries still have non-rejected dependents (rows 1–8, 10, 12–18, 21–28 and 31). | Repoint the deps after re-harvest. Add a check.py warning for any edge into a rejected entry. |
| 47 | 30 `pipeline.*` entries are unreviewed (`pipeline-r1.md` missing; only the `.log` exists). | Rerun the pipeline review. |
| 48 | Symbol-conflict warnings for alpha, lambda and phi point at the rejected `gr.notation-pechenick`. | Drop `conflicts_with` on the rejected side, or have the checker skip rejected entries. |
| 49 | Alias map is missing for these sets: ℏ/hbar; α/alpha_f; σ/σ_SB/σ_sb; R∞/R^inf/R_inf versus the Rydberg R_∞; Φ/Φ(r); Λ/Λ(r); y/y_col/m; g/g_s/g_surf; θ_B (k–B angle) / θ_kB / α_Bk versus surface θ_B (colatitude) versus θ_Bn; E_Bi/E_ci; ω_c/ω_ce; N_H^fit/N_H; P_spin; x_g/u; E_ph/E; B_14 versus B_p,12. Pixel (x,y) also needs scoping against chromaticity (x,y). | One notation-consolidation pass (about 10 entries). |
| 50 | XYZ notation is duplicated in `colour.toml` and `colour_srgb.toml`. | Merge them (with row 27). |
| 51 | 8 blessed notation entries are unreachable from the root: `const.notation-particles`, `const.notation-atomic`, `gr.notation-geometry`, `gr.notation-radiation`, `gr.notation-potekhin`, and `obs.notation-rxj1856-spectral`, `-photometry`, `-timing`. | Wire them to their users or drop them. |
| 52 | The 17 papers fetched today have 0 entries. | They are the sources for 1.3–1.5, 2.3, 2.5 and 3.1–3.9. |
| 53 | Stale rows in refs/INDEX.md. The `_CORRECT` rows name files that no longer exist. The deleted `potekhin_chabrier_2001_eos` row remains. Many entries have venue or volume marked "not in first 60 lines". | Regenerate INDEX.md from CHECKSUMS. |

### 4.4 Provenance "unsure" items (wfA and pass-B), as live decisions

Moot items are already rejected or dropped: gr.visible-mu-exact-values, obs.rxj1856-characteristic-age, obs.approx-vkk2008-dipole-spindown, rt.cond-damping-rates, rt.cond-emissivity-definition, the colour Reinhard white-operator, key-range and zone entries, srgb-gamma22-fit, cmf-origin, opacity.ff-opacity-kramers-gaunt-mcphac, opacity.vacuum-resonance-energy, rt.vacuum-resonance-val, and all eos BSk/TOV items (eos.rho-surface-fe56, m-u-rounded, bsk-adiabatic-index, tov-integration-procedure, the benchmarks).

| # | Entries | Current class | Decision needed |
|---|---------|---------------|-----------------|
| 54 | const.stefan-boltzmann, first- and second-radiation-constant | measured | These are exact given h, k and c. Classify as convention, or derivable and computed (pipeline.constant-set says "computed, not typed"). |
| 55 | const.thomson-cross-section, classical-electron-radius, bohr-radius, the four Rydberg entries, hartree-energy, bohr-magneton | measured | Choose one rule: either the CODATA-adjusted value is measured, or the value is derived from α, m_e and so on. |
| 56 | const.atomic-mass-constant, unified-atomic-mass-unit, vacuum-permeability, vacuum-permittivity | measured | Keep as measured. Note that the first two are definitional. |
| 57 | colour.cie1931-2deg-cmf-file | convention | It is a table, but per A38 it is a definitional standard. Keep convention/leaf and record the reasoning. |
| 58 | colour.srgb-stokes-*-matrix, srgb-inverse-*, srgb-linear-to-xyz-matrix | convention | A39 says derive the matrix. Printed matrices stay as convention verification targets. Add derivable entries for the computed inverses. |
| 59 | colour.reinhard-bounded-output, reinhard-log-average-as-key, reinhard-global-limitation | derivable / convention | Accept as tagged. |
| 60 | gr.beloborodov-cosine-relation, visible-fraction-cosine, cosine-error-bound | fit with drop_reason | Accept: verification targets only. |
| 61 | gr.full-visibility-radius, dphi-max-r16, circle180-radius-r16 | code | Accept as verification targets. Recompute in-house (1.7). |
| 62 | gr.approx-* (schwarzschild-exterior, opaque-surface) | derivable / stub | Approximations need their own convention. Propose `provenance = "convention"` and `ownership = "leaf"` for modelling assumptions. This is a schema decision. |
| 63 | gr.schwarzschild-radius, gr.emission-angle-definition | convention | Accept. |
| 64 | obs.rxj1856-dipole-field | derivable | Reclassify as an inference that depends on assumed I and R. Mark it not used if B is fitted (A35). |
| 65 | obs.rxj1856-extrapolated-stmag-f475w | measured | Reclassify as code or fit: it is an extrapolation of a model. |
| 66 | obs Burwitz fit results (xray-kt, -rinf, -nh, -bb-r-over-d) and the descriptive-model approximations | measured | Fit results to data. Decide between measured and code. Recommend code with reimplement = 2.4–2.6. |
| 67 | obs.rxj1856-ho2007-b-field, obs.rxj1856-ho2007-nh | code | B was set by hand in Ho 2007, not fitted: label it as an assumption. Keep N_H as code. |
| 68 | rt.cond-reflectivity-from-amplitudes, rt.cond-linear-polarisation-degree | convention | Reclassify as derivable/stub. |
| 69 | rt.cond-approx-bc-truncation | convention | This is a silent clamp (rule 5). Require a loud bound in reimplement. |
| 70 | rt.cond-tcrit-fit-fe | fit | Reclassify as code (DFT output). |
| 71 | opacity.field-atomic-b0, opacity.field-qed-bq | derivable | Accept. |
| 72 | opacity.vacuum-resonance-density | fit | Becomes derivable after 3.7 (root find on our own tensor). |
| 73 | opacity.approx-ff-non-born-correction | fit | Correct as fit (Hummer Padé). Reimplement via 3.1. |
| 74 | rt.eddington-temperature-profile-mcphac | fit/stub | Keep. A12 L0 is in 3.8. |
| 75 | rt.eddington-factors-mcphac | convention | Accept. The f = 1/3 and h = 1/2 values are unanchored and need their own entry. |
| 76 | rt.pjump-conversion-spw09, rt.pjump-val | derivable/stub | Accept. The derivation is in 3.7. |
| 77 | surface.gh83-temperature-distribution | derivable | Accept for GH83 eq. 14. The v1 cos² ansatz stays unsourced. |
| 78 | surface.gh83-approx-blackbody-emission, surface.ppp15-fit-flux-recipe | derivable | Accept. |
| 79 | surface.gh83-synthetic-spectrum-fit, surface.gh83-spectrum-conclusion | code | Accept as verification targets. |
| 80 | surface.phc16-rxj1856-polar-field, -polar-ts, -teff-inf | code | Accept: outputs of Ho+2007. Reimplement via A34. |
| 81 | surface.ppp15-approx-fit-validity, surface.ppp15-fit-structure | fit | Accept. Reimplement via 3.9. |
| 82 | eos.relativistic-units-pr-tr | derivable | Accept. |

---

## Counts

| Section | Gaps / rows | Est. entries |
|---------|-------------|--------------|
| 1. Milestone 1 | 13 gaps | ~47 |
| 2. Milestone 2 | 6 gaps | ~29 (+ the deep dust/Mie and Opacity Project routes) |
| 3. Milestone 3 | 12 gaps (one is a decision) | ~122 with the full envelope chain (~99 with a parametrised map), + ~20 if partial ionisation is adopted |
| 4. Hygiene | 82 rows: 35 re-harvest (4.1), 9 disputes (4.2), 9 graph/schema (4.3), 29 provenance decisions (4.4) | ~45 re-proposed/split entries + ~10 notation consolidations |

Next session order:
1. 1.1–1.8, which are almost all local. They need Rybicki & Lightman (a book) for 1.6 and 1.8.
2. Hygiene rows 1, 2, 31, 45 and 47.
3. 1.11 and 1.12.
4. Then M2 and M3 in table order.
