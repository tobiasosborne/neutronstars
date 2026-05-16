# Master List of Approximations

Every approximation used in the pipeline, with name, validity range, source, and impact assessment.

## Phase 2 (Tracer Bullet)

### A1: Modified blackbody atmosphere
**Description:** I_ν = f_col⁴ B_ν(f_col T) cos(θ_e)^p instead of self-consistent RT.
**Source:** Standard hardness-ratio approximation; f_col values from Suleimanov+ (2009).
**Validity:** Good for broadband flux within ~10-20%; poor for spectral features and polarisation.
**Impact:** HIGH — removes all spectral lines, absorption edges, and polarisation dependence. Replaced in Phase 3.
**af node:** 1.3.1

### A2: No rotation (Schwarzschild metric)
**Description:** Ignore NS spin; no Doppler boosting, no frame dragging.
**Source:** Valid when v_rot << c, i.e., P_spin >> R/c ≈ 40 μs.
**Validity:** Good for P > 10 ms (most isolated NSs). Breaks for millisecond pulsars.
**Impact:** LOW for slow rotators. Removes asymmetric brightness due to Doppler boost.
**af node:** 1.5.4

### A3: Fully ionised hydrogen
**Description:** Ignore bound states of H in strong B. All electrons free.
**Source:** Valid when kT >> E_bind(H, B). For B ~ 10¹² G, E_bind ~ 160 eV, so valid for T >> 2×10⁶ K.
**Validity:** Good for T_eff > 5×10⁵ K at B < 10¹³ G (Potekhin & Chabrier 2004).
**Impact:** MEDIUM — at lower T, bound-free opacity creates absorption features. Deferred to Phase 4.

### A4: Local vacuum polarisation without mode conversion
**Description:** Include QED vacuum-polarized local normal-mode vectors in the magnetic opacity weights, but do not yet propagate photons through the vacuum resonance with partial/adiabatic mode conversion.
**Source:** Potekhin, Lai & Chabrier (2004) for vacuum coefficients and mode vectors; van Adelsberg & Lai (2006) for the B-frame cyclic components and atmosphere treatment.
**Validity:** Local opacities away from resonance regions are improved. Spectra across the vacuum resonance are incomplete until mode conversion is implemented.
**Impact:** HIGH for magnetar-strength fields and resonance-region spectra/polarization. The current Fig. 2 metric masks the vacuum-resonance band.

### A5: Plane-parallel atmosphere
**Description:** Each surface element has a plane-parallel atmosphere. Ignore curvature.
**Source:** Standard for NS atmospheres where H_scale << R.
**Validity:** H_scale/R ~ kT/(m_p g R) ~ 10⁻⁷ for typical parameters. Excellent approximation.
**Impact:** NEGLIGIBLE.

### A6: Greenstein-Hartke T(θ) model
**Description:** T(θ_B) = T_eq + (T_pole - T_eq) cos²(θ_B). Parametric, not self-consistent.
**Source:** Greenstein & Hartke (1983).
**Validity:** Qualitatively correct for dipole geometry. Quantitatively approximate — self-consistent thermal transport (Potekhin, Pons & Page 2015) gives different profiles.
**Impact:** MEDIUM — affects hotspot shape and size. Sufficient for tracer bullet.

### A7: Observer at infinity
**Description:** Observer distance D → ∞ so all rays are parallel at the observer.
**Source:** Standard for D >> R (always true for astrophysical observations).
**Validity:** D ~ 100 pc >> R ~ 10 km. Perfect.
**Impact:** NEGLIGIBLE.

## Phase 3+ (to be updated)

### A8: Scalar scattering approximation in magnetic RT
**Description:** Magnetic RT uses local per-mode scattering extinction/albedo, but not the full angle- and mode-redistribution scattering kernel from the published transfer equations.
**Source:** Suleimanov, Potekhin & Werner (2009) Eqs. 5-10 describe mode-coupled scattering source terms.
**Validity:** Useful for flux-normalized smoke tests and local opacity work; not final for physically exact spectra or polarization.
**Impact:** HIGH for oblique-angle and scattering-dominated branches. Likely relevant to the remaining `theta_B=45 deg` Fig. 2 mismatch.

### A9: Ground Landau level only (future)
**Description:** Assume all electrons in ground (n=0) Landau level.
**Validity:** ℏω_c >> kT, i.e., B >> 4.7×10⁸ (T/10⁶)² G.
**Source:** Mészáros (1992) p.47.
**Impact:** Simplifies opacity expressions. Valid for B > 10¹⁰ G at T < 10⁷ K.

### A10: No magnetic pressure contribution to hydrostatic equilibrium (future)
**Description:** B²/(8π) << gas pressure in atmosphere.
**Validity:** B < 10¹⁵ G. For B ~ 10¹² G and T ~ 10⁶ K, P_mag/P_gas ~ 10⁻⁸.
**Source:** To be verified quantitatively.
**Impact:** NEGLIGIBLE for B < 10¹⁴ G.

### A11: True Greenstein-Hartke (1983) T(θ) form — deferred
**Description:** The shipped surface-temperature model (A6) is the
phenomenological cos² ansatz, not the GH 1983 derivation. GH 1983
derive `T(θ) ∝ |cos θ|^{1/4}` for a pure dipole with anisotropic
electron thermal conductivity (their actual Eq. 1 is the dipole *field
magnitude*, not the temperature). A true GH implementation
(`greenstein_hartke_temperature(θ_B, T_pole) = T_pole · |cos θ_B|^{1/4}`,
clamped at a small floor near the equator) is not yet provided.
**Source:** Greenstein & Hartke, ApJ 271, 283 (1983).
**Validity:** GH form is valid for a pure dipole with anisotropic
conductivity in the limit of strong B and where lateral conduction is
negligible. It vanishes at the magnetic equator (unphysical for finite
conductivity).
**Impact:** MEDIUM — would change every rendered hotspot shape. Tier-1
tests assume the cos² form (see A6 / D8). Deferred to a deliberate
physics decision, ideally in tandem with self-consistent thermal
transport (Potekhin, Pons & Page 2015).
**Status:** Not implemented. Tracked here so the GH name doesn't get
re-attached to the cos² ansatz by mistake. See `notes/decisions.md`
D8.
