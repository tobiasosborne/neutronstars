`python3 -I truth/check.py` exited **0**: `318 entries: {'proposed': 318}; 0 errors, 249 warnings`. Warnings: 246 missing provenance classes; 3 GR notation conflicts.

Line references: **H** = [refs/haakonsen_2012_mcphac.txt](/home/tobias/Projects/neutronstars/refs/haakonsen_2012_mcphac.txt), **S** = [refs/suleimanov_potekhin_werner_2009_mag_atm.txt](/home/tobias/Projects/neutronstars/refs/suleimanov_potekhin_werner_2009_mag_atm.txt), **V** = [refs/van_adelsberg_lai_2006_magnetic_atm.txt](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.txt).

| id | verdict | reason |
|---|---|---|
| rt.notation-mcphac | fix | Locator excludes several definitions. Replace with `"Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4"`: H:155–223, 260–293, 359–371, 405–418, 1107–1115. |
| rt.notation-spw09 | fix | Set η units to `"erg s^-1 g^-1 Hz^-1 sr^-1"` and εΛ units to `"erg s^-1 g^-1"` from S:223–230, 463–470. Define φ between projected **B and k**, as Fig. 1 shows; S:111–112 incorrectly says projected n. |
| rt.notation-val | accept | Definitions match V:291–332, 459–478, 517 and 533–534; κT is an opacity, not a cross section. |
| rt.notation-vacuum | fix | Replace EBi meaning with `"ion cyclotron energy"`; remove its unconditional numerical formula. V:224 prints a field-independent value, conflicting with S:267–274, which includes B dependence. Page image confirms the source defect. |
| rt.hydrostatic-column-mcphac | accept | Anchor and equation match H:198–210, including neglected radiation pressure and constant gravity. |
| rt.feautrier-second-order-mcphac | fix | Page image confirms eq. 2 prints denominator `dτν`, without its square (H:151–154). Add eq. 13 to locator and explicitly document the correction: eq. 13 supplies `dτ_j^2` (H:355–360). |
| rt.source-function-mcphac | accept | Thermal coefficient, scattering factor 2, angular integral and albedo match H:169–178, 199–203. |
| rt.feautrier-discrete-mcphac | fix | Replace scattering coefficient `-2\rho_\nu` with `-2\rho_\nu/\sigma_{T,r}`; label this as derived from eqs. 4–6, rather than literal eq. 7. H:174–175 and 1065–1068 establish the missing normalisation; H:213–218 prints the inconsistent form. |
| rt.dipole-phase-function-mcphac | accept | Both kernels match H:299–306 and 324–330; page image confirms 3/16 and every squared term. |
| rt.feautrier-matrix-coeffs-mcphac | reject | Page image confirms eq. 14 prints `1` in the coupling matrix and `dσ/dΩ` (H:355–370). The identity term required by eq. 13 and conversion to the `dμ′` kernel of eq. 4 are unresolved; this cannot be blessed as an implementable matrix. |
| rt.boundary-conditions-mcphac | accept | Conditions and second-order diffusion approximation match H:162–165, 308–317. |
| rt.emergent-intensity-mcphac | fix | Equation matches H:322–323, but set `provenance = "derivable"`: it follows from the Feautrier definition and zero incident intensity, H:158–164. |
| rt.eddington-temperature-profile-mcphac | accept | H:275–293 explicitly prints positive exponent 3 and `0.265 T_eff`; page image agrees. Entry correctly restricts this fit to initialisation and supplies a recomputation route. |
| rt.eddington-factors-mcphac | accept | All three integrals and their normalisations match H:359–367 and 414–418. |
| rt.te-delta-t-mcphac | accept | Numerator, denominator and sign match H:1107–1133; first-order linearisation is stated. |
| rt.flux-constancy-mcphac | accept | Flux definition, 1/(4π), κν weighting and equilibrium zero match H:1081–1106. |
| rt.te-moment-equation-mcphac | fix | Remove eq. 18 from the asserted equivalence and use `"eqs. A6, A26–A27"`. Eq. 18 prints `1/kν` (H:347–357); A6 prints `(1−ρν)` (H:1138–1146). Current latex matches A6 and H:1166–1174. |
| rt.te-rybicki-matrices-mcphac | fix | Set `image_verified = true`: PDF pp. 15–16 confirm A12–A23, H:1136–1154 and 1175–1218. Clarify that `π = Δ1+Δ2` is the paper’s grid variable, not numerical π. |
| rt.observed-flux-redshift-mcphac | fix | Equations match H:226–243. Delete the notes’ `"no light bending"` and `"flat-space solid angle"` claims; replace with `"Uniform-surface distant-observer flux benchmark."` Apparent-radius magnification is explicitly given in refs/potekhin_pons_page_2015_thermal.txt:2106–2111. |
| rt.hydrostatic-spw09 | accept | All gravity, radiative-acceleration and column-depth factors match S:120–145. |
| rt.two-mode-transfer-spw09 | accept | Page image confirms eqs. 5–10: Bν/2, both 1/(2π) factors, indices and angular limits match S:148–230. |
| rt.energy-balance-spw09 | accept | Eq. 9 matches S:211–221, including the outer 1/(2π). |
| rt.surface-flux-spw09 | accept | Eqs. 18–20 match S:367–387, including the direction of the flux-error ratio. Convergence numbers match S:423–425. |
| rt.boundary-conditions-spw09 | accept | Both conditions, including Bν/2, match S:323–333. |
| rt.pjump-conversion-spw09 | accept | Exchange rule, negative exponent, π/2 and cubic power match S:334–365. |
| rt.temperature-correction-spw09 | accept | Page image confirms signs, mode sum, Λ denominator and averaged absorption factor in S:463–487; averaging and layer assignments match S:390–422. Missing Kurucz details are acknowledged. |
| rt.two-mode-transfer-val | fix | Equations match V:281–332. Add validity: `"Normal-mode transfer requires the Faraday depolarisation condition; it breaks down near the vacuum resonance."` This limitation is explicit at V:268–274 and 585–588. |
| rt.flux-constraints-val | accept | Both constraints match V:297–332, including Bν/2, 2π and the outgoing-minus-incoming sign. |
| rt.source-function-val | reject | Fresh extraction and the p. 7 image confirm that eq. 30 really contains mixed mode indices and unprimed μ inside the μ′ integral (V:484–493). Literal transcription is confirmed, but the general scattering source function is defective; eq. 32 does not repair that bundled claim. |
| rt.boundary-conditions-val | accept | Limits, Bν/2, zero incoming intensity and 5–8 decades match V:541–547. |
| rt.vacuum-resonance-val | fix | Add validity `"E \ll E_{Be}"`, the regime introduced at V:173–187. Coefficients, exponents and units otherwise match V:181–186 and 213–215. |
| rt.adiabatic-energy-val | fix | Eq. 3 matches the p. 4 image and V:216–225. Replace the notes’ EBi numerical assertion with `"The printed E_Bi definition omits field dependence; use a separately sourced cyclotron-energy definition."` Compare S:267–274. |
| rt.pjump-val | accept | Probability, both exchanges and limiting cases match V:227–232, 251–255 and 529–545. |
| rt.vacuum-relevance-val | accept | Scaling powers and the distinction between spectral and polarisation effects match V:235–244. |
| rt.no-feautrier-with-conversion-val | fix | Rename to `"Standard Feautrier substitution not used for VAL06 partial conversion"` and restrict validity to that procedure. Remove the low-B exemption: V:536–545 describes this implementation limitation; V:242–244 retains polarisation effects below Bl. |
| rt.mode-stokes-relations-val | fix | Add `"neglecting dissipative terms in the dielectric tensor"` to the statement, as required by V:615–622. All signs and factors in eqs. 43–47 match V:624–638. |
| rt.approx-plane-parallel-mcphac | accept | Scale-height and constant-gravity justifications match H:183–185 and 203–210. |
| rt.approx-coherent-scattering-mcphac | fix | Replace validity with `"Coherent-scattering assumption; inaccuracies may occur above 5 keV and at very high T_eff; no guaranteed safe cutoff is supplied."` H:174–177 and 1005–1009 give warnings, not a validated below-5-keV domain. |
| rt.approx-no-radiation-pressure-mcphac | accept | H:186–187 explicitly gives `T_eff ≪ 10^7 K`; p. 3 image confirms the inequality obscured in extraction. |

- Canonical notation remains unresolved: RT μ is a positive half-range variable in H/V, signed in S, and differs from GR μ = cos ψ. Photon ν also collides with timing ν; McPHAC j labels angular points while magnetic-transfer j labels modes.
- Merge documented aliases explicitly: y/y_col/m; g/g_s/g_surf; θkB/αBk/opacity θB; EBi/Eci; and the duplicated vacuum symbols. Greek, ASCII and LaTeX spellings currently conceal matches.
- V:186 delegates the actual function fB to Ho & Lai (2003), absent locally. These entries establish neither a computable fB nor permission to substitute fB = 1.
- All 39 entries have a plausible Milestone 1–3 role; no clear `drop`.