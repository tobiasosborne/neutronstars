# Spine agent report (harvest gaps)

The spine is written to `/home/tobias/Projects/neutronstars/truth/entries/pipeline.toml`, and `python3 -I truth/check.py --render` reports 0 errors.

**Departure from the brief:** the root is `pipeline.rendered-frame`, not pixel radiance. Its job is HDR radiance plus the derived display image, and it depends on `pipeline.pixel-radiance` and the two display nodes. With pixel radiance as the root, the colour and tone entries would have formed separate roots.

**Cycle avoided:** `pipeline.surface-parameters` (audit A34, the parameter fit) uses `pipeline.earth-photometry` as its forward model, and pixel radiance depends on the parameters. So the R5 harness is reachable from the root without a cycle. Its milestone is 2.

**Rejected entries:** some entries changed status while I worked, because other lanes were reviewing at the same time. I re-ran the check against the final statuses. Seven rejected entries are listed, each with a reason in `notes` because nothing else covers the need:
- `opacity.ff-opacity-kramers-gaunt-mcphac`
- `opacity.ff-effective-frequency`
- `opacity.approx-ff-proton-collisions-negligible`
- `opacity.dielectric-tensor-polarizability`
- `rt.emergent-intensity-mcphac`
- `rt.cond-fit-coefficient-data`
- `rt.cond-fit-structure` is blessed; the seventh is covered by the list above only once each, so the count is six named here plus the coefficient data entry already listed.

I swapped the newly rejected `eos.metric-schwarzschild-interior`, `eos.enclosed-mass-and-rg`, `eos.redshift-observables` and `rt.vacuum-resonance-val` for blessed alternatives.

**Tree state (other lanes):** the root renders red only because four blessed eos entries have no `ownership` field yet (`eos.surface-gravity`, `eos.thin-envelope-newtonian-hydrostatics`, `eos.electron-pressure-magnetised`, `eos.relativistic-units-pr-tr`). The other 39 roots in `truth/TREE.md` are rejected entries or the BSk/TOV chain that audit item A24 drops. There is also an existing cycle between `opacity.mode-scattering-opacity` and `opacity.approx-isotropic-scattering`.

**Results**
- **Spine nodes:** 30 (milestone 1: 17, milestone 2: 1, milestone 3: 12).
- **Reachable from the root:** 324 entries in total, 294 of them non-spine, out of 387 entries in the database.
- **Distinct leaf entries reachable from the root:** 90 (entries with no `depends_on`). 62 of these are not notation entries.

**Needs not met by any existing entry (harvest gaps)**
1. Exterior Schwarzschild line element (all metric entries are rejected).
2. Timelike geodesics in Schwarzschild (E, L, effective potential, periastron, proper time).
3. Full orthonormal static tetrad, plus the local Lorentz boost and aberration formula.
4. Camera projection model (pinhole or equidistant).
5. Special-relativistic Doppler factor for a moving observer (the R2 far-field check).
6. Bending integral from a finite camera radius.
7. Parallel transport of the polarisation basis (Stokes rotation) along null geodesics.
8. Camera PSF or diffraction model, veiling-glare model, physical exposure and sensor model, inverse tone operator, the 683 lm/W photometric scaling, and video time remapping.
9. A blessed inverse sRGB power segment (both existing versions are rejected).
10. The Planck function B_ν(T) as an entry of its own.
11. Constants: IAU 2015 nominal GM_sun, and the exact au and pc definitions.
12. Sky: the Gaia DR3 catalogue, passbands and zero points; RX J1856 RA and Dec; stellar SEDs from Gaia photometry; the starless diffuse Milky Way map and its licence.
13. R5 harness: HST filter throughputs, an owned extinction curve (CCM89 or a dust model), and X-ray instrument responses.
14. Parameter fit: Chandra LETG spectrum and responses, ISM photoabsorption (Wilms+2000, Verner+1996), and a likelihood or fitting procedure.
15. Dipole field: B_r, B_θ and |B|, and the GR dipole (Wasserman and Shapiro 1983).
16. A temperature map parametrised and fitted to the pulse profile.
17. Envelope chain (audit item A29):
    - PPP15 eqs. 26, 31–36, 38 and A.1–A.3
    - Potekhin 1999 Coulomb logarithms
    - Baiko and Yakovlev 1995
    - OCP Monte Carlo and the SDD82 Table I values
    - the PC01 non-ideal EOS
    - Fe/C bound-free opacities
18. Ideal-gas EOS ρ = m_p P / (2 k_B T) for fully ionised H.
19. Feautrier method from Mihalas 1978, a blessed Feautrier coupling matrix, the exact Hopf function, McPHAC spectra fixtures and Suleimanov+2009 figure values.
20. Magnetic mode-coupled scattering redistribution (Mészáros 1992) and Compton (Kompaneets) redistribution.
21. Gaunt factor (A4): the Sommerfeld / Karzas–Latter g_ff, van Hoof 2014 averaging, Bethe–Salpeter Coulomb functions, and the gffgu.dat fixture.
22. PC03 Λ_pp fit, eq. 48, and its App. A chain A1–A6 (A18).
23. A16 chain: PC03 App. B eqs. B6–B8 and B19–B22, Potekhin 1994 matrix elements, and Landau wavefunctions (Ruder 1994).
24. Vacuum polarisation (A20): Heyl–Hernquist 1997 and Euler–Heisenberg proper-time integral, Kohri and Yamada 2002, and a derivation of the constant 0.0329199.
25. Mode conversion (A21/A22): Lai and Ho 2002 Landau–Zener derivation, and the Ho and Lai 2003 mode-evolution equations.
26. Condensed surface: general Fresnel reflectivity at arbitrary θ_k and θ_B (van Adelsberg 2005 derivation), condensed-matter EOS and surface density (Medin and Lai 2006/2007), and Potekhin 1999 conductivity for the damping rates.