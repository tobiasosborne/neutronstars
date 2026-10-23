`python3 -I truth/check.py` exited **0**:
```text
357 entries: {'blessed': 220, 'rejected': 46, 'proposed': 91}; 0 errors, 524 warnings
```
Warnings: 208 missing provenance, 313 missing ownership, 3 GR notation conflicts.

Line references: [H = McPHAC](/home/tobias/Projects/neutronstars/refs/haakonsen_2012_mcphac.txt), [S = SPW09](/home/tobias/Projects/neutronstars/refs/suleimanov_potekhin_werner_2009_mag_atm.txt), [V = VAL06](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.txt), [L = Ho & Lai 2003](/home/tobias/Projects/neutronstars/refs/ho_lai_2003_atmospheres_vacuum.txt).

| id | verdict | reason |
|---|---|---|
| rt.notation-mcphac | accept | Anchor matches H:160. Expanded locator covers the definitions: H:155–223, 260–293, 359–371, 414–418, 1107–1115. |
| rt.notation-spw09 | accept | Anchor matches S:138. Corrected emissivity and energy-error units follow S:223–230, 463–470. PDF Fig. 1 confirms φ between projected B and k; S:111–112 contains the acknowledged prose typo. |
| rt.notation-vacuum | accept | Anchor matches V:224. Remaining definitions match V:181–247; the defective numerical ion-cyclotron assertion has been removed. |
| rt.feautrier-second-order-mcphac | accept | Anchor matches H:152. PDF images confirm the missing square in printed Eq. 2 and the squared denominator in Eq. 13, H:355–358. Correction and provenance are explicit. |
| rt.feautrier-discrete-mcphac | accept | Anchor matches H:203. The corrected coefficient \(2\rho_\nu/\sigma_{T,r}=2/(\kappa_\nu+\sigma_{T,r})\) follows H:169–203 and agrees with A1, H:1065–1068. It is explicitly distinguished from defective printed Eq. 7. |
| rt.emergent-intensity-mcphac | reject | Locator incorrectly says “after eq. 17 discussion”; the equation **precedes** Eq. 17, H:322–330, confirmed visually. Replace locator with `"p. 5, Sect. 3.3, unnumbered equation before eq. 17"`. **Surviving weaker statement:** \(I_\nu(0,\mu)=2P_\nu(0,\mu)\) for zero incident intensity, H:162–164, 322–323; `derivable` is supported. |
| rt.te-moment-equation-mcphac | accept | Anchor matches H:1147. A6, including the entire bracket multiplied by \(1-\rho_\nu\), matches H:1134–1155; both boundary conditions match H:1166–1174. PDF image confirms grouping. |
| rt.te-rybicki-matrices-mcphac | accept | Anchor matches H:1216. Fresh extraction and visual inspection of PDF pp. 15–16 confirm A12–A23, including factors 8, negative coupling, and grid variable \(\pi=\Delta_1+\Delta_2\): H:1136–1154, 1175–1218. |
| rt.observed-flux-redshift-mcphac | accept | Anchor and both equations match H:226–243: inverse first power of \(1+z\), \(R^2/D^2\), and blueshifted frequency argument. Attribution to McClintock et al. is retained. |
| rt.two-mode-transfer-val | accept | Anchor and equations match V:281–332, including directional signs and \(\kappa_j^{tot}/\kappa_T\). Added Faraday limitation matches V:268–274, 585–588. |
| rt.vacuum-resonance-val | reject | Anchor, equations and added regime match V:173–187, 213–215, but the current note claiming Ho & Lai (2003) is “not in refs/” is false: L:444–454 supplies Eq. 2.41 locally. Delete that assertion and link the local definition. **Surviving weaker statement:** Eqs. 1–2 hold in the stated \(E\ll E_{Be}\) regime with the referenced \(f_B\); the missing-source claim does not survive. |
| rt.adiabatic-energy-val | accept | Anchor matches V:216. PDF p. 4 confirms absolute-value bars, exponents \(2/3\) and \(1/3\), and coefficient 2.52, V:216–225. Implicit energy dependence is acknowledged, V:247. |
| rt.no-feautrier-with-conversion-val | accept | Anchor matches V:536. Restricted claim and resonance-stop integration procedure faithfully reproduce V:536–545; no general impossibility or low-field exemption is asserted. |
| rt.mode-stokes-relations-val | accept | Anchor matches V:627. Dissipation qualification matches V:615–622; every sign and factor in Eqs. 43–47 matches V:624–638. |
| rt.approx-coherent-scattering-mcphac | accept | Anchor and energy warning match H:174–177; frequency decoupling matches H:317–326. High-temperature caveat matches H:1005–1009. No validated safe cutoff is claimed. |

- Scope μ explicitly: H:155–160 and V:291 use a positive half-range variable; S:130–140, 148–209 uses signed directions. McPHAC angular index \(j\) also differs from magnetic mode index \(j\); frequency index \(k\) differs from opacity \(k_\nu\).
- Consolidate remaining aliases: \(y/y_{\rm col}/m\), \(g/g_s/g_{\rm surf}\), \(\alpha_{Bk}/\theta_{kB}/\text{opacity }\theta_B\), and Greek/ASCII vacuum-symbol spellings.
- Wire \(E_{Bi}\) to a separately sourced cyclotron definition: V:224 omits field dependence; S:267–274 supplies it. All 15 entries have a Milestone 2–3 role.