# Notation dictionary (rendered by truth/check.py; do not edit)

| symbol | meaning | units | entry | source |
|---|---|---|---|---|
| \alpha | emission angle: angle between the photon direction and the local radial direction (surface normal), measured by the local static observer at the emission point | rad | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| \alpha_f | fine-structure constant e^2/(hbar c) | 1 | const.notation-atomic | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45–025010-46) |
| \alpha_{opt} | (paper symbol alpha) optical/UV power-law index, F_lambda ~ lambda^-alpha (Rayleigh-Jeans: 4) | 1 | obs.notation-rxj1856-photometry | kaplan_2011_rxj1856_spectrum (Sec. 3.1 and Table 5; Sec. 2.3 and Table 2) |
| \alpha_{RA} | right ascension of a source; the reference frame and epoch are attributes of each datum (Gaia DR3: ICRS, Gaia-CRF3, epoch J2016.0; Walter+2010: not stated) | rad | obs.notation-gaia-sky | lindegren_2021_gaia_astrometry (Sect. 3.1, p. 6 (refs/lindegren_2021_gaia_astrometry.txt:432-437,470-475)) |
| \bar{x} | CIE 1931 2-deg colour matching function x-bar | source-unspecified (tabulated values) | colour.notation-cie-cmf | data/cie1931_2deg_cvrl_README (line 2 (Data: ... columns)) |
| \bar{y} | CIE 1931 2-deg colour matching function y-bar | source-unspecified (tabulated values) | colour.notation-cie-cmf | data/cie1931_2deg_cvrl_README (line 2 (Data: ... columns)) |
| \bar{z} | CIE 1931 2-deg colour matching function z-bar | source-unspecified (tabulated values) | colour.notation-cie-cmf | data/cie1931_2deg_cvrl_README (line 2 (Data: ... columns)) |
| \beta | bending angle beta = psi - alpha | rad | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| \chi_n | chi_n = (mu - m_e c^2 sqrt(1+2bn))/(k_B T), reduced chemical potential of level n | 1 | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| \cos\psi | mu = cos psi | 1 | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| \cos\psi_{vis} | critical value of mu = cos psi at which the emission is tangential (alpha = 90 deg); surface points with mu < mu_v are invisible (single-image regime R > 2 r_g only; compact stars with R < 2 r_g or R < 3 GM/c^2 show multiple images) | 1 | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| \Delta T | temperature correction Delta T(y), T(y) = T_0(y) + Delta T(y) | K | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| \Delta T_\Lambda | temperature correction from the integral Λ-iteration | K | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \Delta\theta_{lt} | angular-position error caused by neglecting light-time effects in the epoch propagation (first-order estimate) | rad | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \delta_{Dec} | declination of a source; the reference frame and epoch are attributes of each datum (Gaia DR3: ICRS, Gaia-CRF3, epoch J2016.0; Walter+2010: not stated) | rad | obs.notation-gaia-sky | lindegren_2021_gaia_astrometry (Sect. 3.1, p. 6 (refs/lindegren_2021_gaia_astrometry.txt:432-437,470-475)) |
| \delta_{L} | small offset added to luminance to avoid the singularity at black pixels (paper symbol delta) | same as L_w | colour.notation-reinhard-parameters | reinhard_2002_tonemapping (p. 4, eq. 1-2 text) |
| \dot P | period derivative | s s^-1 | obs.notation-rxj1856-timing | vankerkwijk_kaplan_2008_rxj1856_timing (Table 2 note) |
| \dot\nu_{spin} | spin frequency derivative (paper symbol nu-dot) | Hz s^-1 | obs.notation-rxj1856-timing | vankerkwijk_kaplan_2008_rxj1856_timing (Table 2 note) |
| \eta_\nu^j | emissivity in mode i, per unit mass-column normalisation of eq. 10 | erg s^-1 g^-1 Hz^-1 sr^-1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \eta_{\rho} | dimensionless factor eta in the surface density rho_s: uniform-electron (Wigner-Seitz) model, unity for the ion-sphere model | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| \Gamma_{coul} | Coulomb coupling parameter Gamma = (Ze)^2/(a_ion k_B T) of the ions (van Adelsberg+2005 sect. 2.1) | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| \gamma_{ei}^{(l)} | dimensionless collisional damping rate hbar nu/E for motion longitudinal to B | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| \gamma_{ei}^{(tr)} | dimensionless collisional damping rate hbar nu/E for motion transverse to B | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| \kappa^a_j | true-absorption (free-free for fully ionised H) opacity of normal mode j; depends on direction through theta_kB | cm^2 g^-1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \kappa^a_j | true-absorption (free-free for fully ionised H) opacity of normal mode j; depends on direction through theta_kB | cm^2 g^-1 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| \kappa^s_j | total scattering opacity out of normal mode j, kappa^s_{j1} + kappa^s_{j2}; depends on direction through theta_kB | cm^2 g^-1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \kappa^s_j | total scattering opacity out of normal mode j, kappa^s_{j1} + kappa^s_{j2}; depends on direction through theta_kB | cm^2 g^-1 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| \kappa^s_{jj'} | scattering opacity kappa^s_{jj'} from mode j into mode j' | cm^2/g | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| \kappa_\nu | B = 0 true-absorption opacity at photon frequency nu (free-free for fully ionised hydrogen; McPHAC also allows other absorption processes) | cm^2 g^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| \kappa_j | total extinction opacity of normal mode j, kappa^a_j + kappa^s_j | cm^2 g^-1 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| \kappa_T | Thomson scattering opacity, 0.4 cm^2 g^-1 | cm^2 g^-1 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| \lambda | wavelength column of the CMF file | nm | colour.notation-cie-cmf | data/cie1931_2deg_cvrl_README (line 2 (Data: ... columns)) |
| \Lambda | metric function Lambda(r): g_rr = e^(2 Lambda) | 1 | gr.notation-potekhin | potekhin_pons_page_2015_thermal (Sect. 2.1 (eqs. 1-3) and Sect. 5 (eqs. 39-40)) |
| \Lambda^j_{\nu,diag} | diagonal element of the Λ operator for mode i, a function of the mean optical depth | 1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \lambda_e | electron Compton wavelength hbar/(m_e c) | cm | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| \lambda_{aff} | affine parameter along the photon orbit (normalisation source-specific: Bel02 lambda -> r as r -> infinity; Pech83 eq. 2.2) | cm | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| \langle F_\lambda\rangle | synthetic mean energy flux per wavelength unit through a passband | erg s^-1 cm^-2 cm^-1 | obs.notation-gaia-photometry | riello_2021_gaia_photometry (Sect. 7, Eq. (7), p. 15 (refs/riello_2021_gaia_photometry.txt:931-940); Eqs. (13)-(17), p. 17; Table 3, p. 18) |
| \langle F_\nu\rangle | synthetic mean energy flux per frequency unit through a passband | erg s^-1 cm^-2 Hz^-1 | obs.notation-gaia-photometry | riello_2021_gaia_photometry (Sect. 7, Eq. (7), p. 15 (refs/riello_2021_gaia_photometry.txt:931-940); Eqs. (13)-(17), p. 17; Table 3, p. 18) |
| \mathrm{au} | astronomical unit, exactly 149 597 870 700 m (IAU 2012 B2) | cm | obs.notation-gaia-sky | lindegren_2021_gaia_astrometry (Sect. 3.1, p. 6 (refs/lindegren_2021_gaia_astrometry.txt:432-437,470-475)) |
| \mathrm{KDC} | black digital count | 1 (counts) | colour.notation-srgb-digital | stokes_1996_srgb (eq. 1.3-1.6) |
| \mathrm{Ry} | Rydberg energy | J | const.notation-atomic | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45–025010-46) |
| \mathrm{STMAG} | magnitude in the STMAG system | mag | obs.notation-rxj1856-photometry | kaplan_2011_rxj1856_spectrum (Sec. 3.1 and Table 5; Sec. 2.3 and Table 2) |
| \mathrm{WDC} | white digital count | 1 (counts) | colour.notation-srgb-digital | stokes_1996_srgb (eq. 1.3-1.6) |
| \mu | cosine of the angle between a ray and the outward surface normal (positive); ingoing rays are written with -μ | 1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| \mu_\delta | proper motion in declination at J2016.0 (catalogue field pmdec) | rad s^-1 | obs.notation-gaia-sky | lindegren_2021_gaia_astrometry (Sect. 3.1, p. 6 (refs/lindegren_2021_gaia_astrometry.txt:432-437,470-475)) |
| \mu_B | Bohr magneton: SI e hbar/(2 m_e); Gaussian e hbar/(2 m_e c) = E_ce/(2B) | J T^-1 | const.notation-atomic | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45–025010-46) |
| \mu_e | electron chemical potential including rest energy | erg | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| \mu_q | Gauss-Legendre angle nodes on (0,1); M of them | 1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| \mu_{\alpha*} | proper motion in right ascension times cos(delta), at J2016.0 (catalogue field pmra) | rad s^-1 | obs.notation-gaia-sky | lindegren_2021_gaia_astrometry (Sect. 3.1, p. 6 (refs/lindegren_2021_gaia_astrometry.txt:432-437,470-475)) |
| \mu_{r} | radial proper motion mu_r = v_r varpi / A (radial velocity in the unit of proper motion), sixth astrometric parameter | rad s^-1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \mu_{tot,0} | total tangential proper motion at the reference epoch, mu_0 = \|mu_0\| = (mu_alpha*0^2 + mu_delta0^2)^(1/2) of BL; spelled mu_{tot,0} here because bare mu_0 is the vacuum permeability (const.notation-fundamental) and mu is the ray cosine | rad s^-1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \mu_{tot} | total tangential proper motion, (mu_alpha*^2 + mu_delta^2)^(1/2) = v_t / b | rad s^-1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \nu | photon frequency in the local rest frame of the emitting plasma (static frame at the surface for a non-rotating star) | Hz | gr.notation-radiation | beloborodov_2002_ray_tracing (Appendix, eq. 11 and text above it) |
| \nu^\infty | photon frequency measured at the observer (at infinity) | Hz | gr.notation-radiation | beloborodov_2002_ray_tracing (Appendix, eq. 11 and text above it) |
| \nu_{ei}^{(l)} | electron-ion collision frequency longitudinal to B | s^-1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| \nu_{ei}^{(tr)} | electron-ion collision frequency transverse to B | s^-1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| \nu_{spin} | spin frequency (paper symbol nu) | Hz | obs.notation-rxj1856-timing | vankerkwijk_kaplan_2008_rxj1856_timing (Table 2 note) |
| \omega_{ce} | electron cyclotron angular frequency omega_ce = eB/(m_e c) (Gaussian CGS) | s^-1 | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| \phi | azimuthal angle of a surface element about the observer direction D | rad | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| \Phi | metric function Phi(r): g_tt = -e^(2 Phi) c^2 | 1 | gr.notation-potekhin | potekhin_pons_page_2015_thermal (Sect. 2.1 (eqs. 1-3) and Sect. 5 (eqs. 39-40)) |
| \Phi_B | angle between the magnetic field lines and the local horizontal (surface tangent plane) (GH83 Fig. 1, Phi) | rad | surface.notation-geometry | greenstein_hartke_1983_surface_T (Figs. 1, 3; eqs. 15, 21) |
| \phi_{rot} | pulsar rotational phase, phi = Omega t (GH83 eq. 21 and Fig. 3 caption) | rad | surface.notation-geometry | greenstein_hartke_1983_surface_T (Figs. 1, 3; eqs. 15, 21) |
| \psi | polar angle of the photon position in its orbit plane, measured from the escape direction (psi = 0 at infinity); for a surface point, the polar angle of the emission point measured from the direction to the observer | rad | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| \Psi_{pole} | geographic colatitude (angle from the rotation axis) of a magnetic pole (GH83 Fig. 3, Psi) | rad | surface.notation-geometry | greenstein_hartke_1983_surface_T (Figs. 1, 3; eqs. 15, 21) |
| \psi_{pole} | angle between the normal to a magnetic pole and the line of sight to the observer (GH83 Sect. III a; written alpha) | rad | surface.notation-geometry | greenstein_hartke_1983_surface_T (Figs. 1, 3; eqs. 15, 21) |
| \rho | mass density (TOV/EOS: gravitational mass density = energy density/c^2; atmosphere/opacity: plasma mass density, equal to it up to O(k_B T/(m_u c^2))) | g cm^-3 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \rho_1 | rho/(1 g cm^-3), dimensionless density | 1 | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| \rho_V | density at which a photon of energy E meets the vacuum resonance | g cm^-3 | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| \rho_{s,cond} | zero-pressure density of the condensed surface (paper: rho_s) | g cm^-3 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| \sigma_T | Thomson cross section (8 pi/3) r_e^2, r_e = e^2/(m_e c^2) | m^2 | const.notation-atomic | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45–025010-46) |
| \sigma_{SB} | Stefan-Boltzmann constant | erg cm^-2 s^-1 K^-4 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| \sigma_{SB} | Stefan-Boltzmann constant | erg cm^-2 s^-1 K^-4 | surface.notation-observer-frame | potekhin_ho_chabrier_2016_review (Sect. 2.1, p. 4; Sect. 2.3, p. 8) |
| \sigma_{T,r} | reduced Thomson scattering opacity of McPHAC (per gram, includes the ionisation fraction) | cm^2 g^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| \tau_\nu | effective normal optical depth, dτ_ν = k_ν dy | 1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| \tau_\nu^j | optical depth in mode i, dτ_ν^i = (k_ν^i+σ_ν^i) dm | 1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \tau_n | tau_n = (1/sqrt(1+2bn)) T/T_r, reduced temperature of level n | 1 | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| \tau_T | Thomson optical depth, dτ = -ρ κ_T dz | 1 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| \tau_{A} | light-travel time of one astronomical unit, tau_A = A/c = 499.004 784 s | s | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \tau_{char} | characteristic age P/(2 Pdot) | Myr | obs.notation-rxj1856-timing | vankerkwijk_kaplan_2008_rxj1856_timing (Table 2 note) |
| \tau_{lt} | light-travel time T_obs - T_em from the source to the barycentre; tau_0 = b_0/c at the reference epoch | s | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \theta | polar angle between the photon wave vector k and the local outward surface normal n; mu = cos theta (Pot12 theta_k is a folded incidence parameter: it equals theta for outgoing rays and pi - theta for incident rays) | rad | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \theta_B | magnetic colatitude: polar angle on the stellar surface measured from the magnetic axis (written beta in GH83 eq. 15 and theta in PPP15 App. B) | rad | surface.notation-geometry | greenstein_hartke_1983_surface_T (Figs. 1, 3; eqs. 15, 21) |
| \theta_{Bn} | angle between the magnetic field B and the local surface normal n | deg | obs.notation-rxj1856-spectral | ho_2007_rxj1856_thin_atm (Sec. 1; Sec. 2; Tables 1-2; Sec. 5.1) |
| \theta_{Bn} | angle between the magnetic field B and the local surface normal n | rad | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \theta_{kB} | angle between the photon wave vector k and the magnetic field B | rad | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \theta_{kB} | angle between the photon wave vector k and the magnetic field B | rad | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| \theta_{los} | colatitude of the line of sight to the observer measured from the rotation axis (written theta in GH83 Fig. 3) | rad | surface.notation-geometry | greenstein_hartke_1983_surface_T (Figs. 1, 3; eqs. 15, 21) |
| \theta_{Sch} | Schwarzschild polar coordinate (orbit plane taken as theta = 90 deg) | rad | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| \varepsilon_\Lambda | energy-balance error as a function of depth | erg s^-1 g^-1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \varepsilon_\parallel | longitudinal element (along B) of the dielectric tensor I + 4 pi chi, frame z \|\| B: eps_zz = eps_parallel | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| \varepsilon_\perp | transverse diagonal element of the dielectric tensor I + 4 pi chi, frame z \|\| B: eps_xx = eps_yy = eps_perp | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| \varepsilon_\wedge | off-diagonal (gyrotropic) element of the dielectric tensor I + 4 pi chi, frame z \|\| B, with eps_xy = +i eps_wedge and eps_yx = -i eps_wedge (convention common to vAL05 and PLCH04) | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| \varepsilon_H | relative flux error 1 - H_0/∫H_ν dν | 1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \varepsilon_{T} | small parameter f_T - 1 of the light-time expansion | 1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \varepsilon_{V} | small parameter f_V - 1 of the light-time expansion | 1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \varphi_k | azimuthal angle of the photon wave vector about the surface normal, measured from the plane containing B and the normal | rad | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \varpi_\nu | dimensionless scattering albedo σ_T,r/(κ_ν+σ_T,r) = σ_T,r/k_ν (not mass density) | 1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| \varpi_\nu^j | scattering albedo of mode i from averaged opacities, σ̄/(k̄+σ̄) | 1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| \varpi_{plx} | trigonometric parallax of a source at the reference epoch J2016.0 (Gaia astrometric parameter) | rad | obs.notation-gaia-sky | lindegren_2021_gaia_astrometry (Sect. 3.1, p. 6 (refs/lindegren_2021_gaia_astrometry.txt:432-437,470-475)) |
| \vec{\mu}_{pm} | proper motion vector du/dt of the barycentric direction | rad s^-1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \vec{b}_{bary} | barycentric position vector of the source: the apparent position b_app(T_obs) = b_true(T_em) from Sect. 4 onward of Butkevich and Lindegren (2014), where the subscripts are omitted; the true position at the same instant is a different quantity, and Sects. 1-3 write b for the true vector where the uniform rectilinear model is stated | cm | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \vec{p}_{tri} | unit vector in the direction of increasing right ascension (first axis of the normal triad) | 1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \vec{q}_{tri} | unit vector in the direction of increasing declination (second axis of the normal triad) | 1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \vec{r}_{tri} | fixed unit vector coinciding with the barycentric direction at the reference epoch; third axis of the normal triad | 1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \vec{u}_{bary} | barycentric coordinate direction to the source, unit vector b/\|b\| in the ICRS | 1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \vec{v}_{app} | apparent space velocity expressed through the astrometric parameters, v = (A_V/varpi)(r mu_r + mu), BL eq. (30); distinct from the true velocity \vec{v}_{true} | cm s^-1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \vec{v}_{bary} | space velocity of the source relative to the solar system barycentre: the constant TRUE velocity in the uniform rectilinear model (BL Sect. 3, eq. 1); the APPARENT velocity v_app from Sect. 4 onward, where the subscripts are omitted. True and apparent velocities differ by the light-time factor 1/(1 - v_r^app/c) of BL eq. (19); never mix them | cm s^-1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| \vec{v}_{true} | true (constant) space velocity of the uniform rectilinear model, BL eq. (1); distinct from the apparent velocity \vec{v}_{app} | cm s^-1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| a_0 | Bohr radius | m | const.notation-atomic | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45–025010-46) |
| a_m | magnetic length (hbar c/(eB))^(1/2) = lambda_e / sqrt(b) | cm | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| a_q | Gauss-Legendre weights on (0,1); M of them | 1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| A_V | visual extinction | mag | obs.notation-rxj1856-photometry | kaplan_2011_rxj1856_spectrum (Sec. 3.1 and Table 5; Sec. 2.3 and Table 2) |
| A_{ion} | ion mass number A | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| a_{ion} | ion-sphere radius a_i, defined by (4 pi a_i^3/3)^(-1) = n_i (ion number density) | cm | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| a_{key} | key value, key value; target scaled luminance in eq. 2, before tone compression; paper symbol a; default 0.18 | 1 | colour.notation-reinhard-parameters | reinhard_2002_tonemapping (p. 4, eq. 1-2 text) |
| A_{vel} | astronomical unit expressed in km yr s^-1: A_V = 4.740 470 446 | km yr s^-1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| b | impact parameter of the photon orbit | cm | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| B'_{sRGB} | non-linear (encoded) sRGB blue value | 1 | colour.notation-srgb-nonlinear | stokes_1996_srgb (eq. 1.2a, 1.2b) |
| B_\nu(T) | blackbody (Planck) intensity at temperature T | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| B_l | approximate, energy- and temperature-dependent field scale above which the vacuum resonance lies between the X- and O-mode photospheres for typical photons (not a universal field threshold) | G | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| B_{12} | magnetic field B/(10^12 G) | 1 | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| B_{14} | local magnetic field strength B/(10^14 G) | 1 | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| B_{8bit} | 8-bit digital code value of the blue channel | 1 (counts) | colour.notation-srgb-digital | stokes_1996_srgb (eq. 1.3-1.6) |
| b_{bary} | barycentric distance of the source, varpi = A/b | cm | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| B_{dip} | dipole magnetic field inferred from spin-down | G | obs.notation-rxj1856-timing | vankerkwijk_kaplan_2008_rxj1856_timing (Table 2 note) |
| b_{rel} | relativistic field parameter b = B/B_Q = hbar omega_ce/(m_e c^2) = B_12/44.14 | 1 | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| B_{sRGB} | linear sRGB blue tristimulus value | 1 | colour.notation-srgb-linear | stokes_1996_srgb (eq. 1.1) |
| cR_{Ry} | Rydberg frequency | Hz | const.notation-atomic | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45–025010-46) |
| D | distance from the observer to the star (D >> R) | cm | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| D | distance from the observer to the star (D >> R) | pc | obs.notation-rxj1856-spectral | ho_2007_rxj1856_thin_atm (Sec. 1; Sec. 2; Tables 1-2; Sec. 5.1) |
| d\kappa^s_{ji}/d\Omega | differential scattering opacity from mode j direction (mu',phi') into mode i direction (mu,phi); under SPW09 definitions dkappa^s_{ji}/dOmega = sigma_nu^{ij}/(2 pi), not sigma_nu^{ij} itself | cm^2 g^-1 sr^-1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| d\Omega | solid angle subtended by dS on the observer sky, dOmega = b db dphi / D^2 | sr | gr.notation-radiation | beloborodov_2002_ray_tracing (Appendix, eq. 11 and text above it) |
| d\sigma_{T,r}/d\mu'(\mu,\mu') | differential amount of scattering from directions μ′ into directions μ | cm^2 g^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| dF | observed flux from the surface element dS, dF = I dOmega | erg s^-1 cm^-2 (frequency-integrated) | gr.notation-radiation | beloborodov_2002_ray_tracing (Appendix, eq. 11 and text above it) |
| dS | proper surface element of the stellar surface, dS = R^2 d(cos psi) d phi | cm^2 | gr.notation-radiation | beloborodov_2002_ray_tracing (Appendix, eq. 11 and text above it) |
| E | photon energy hbar omega | keV | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| E_1 | photon energy E / (1 keV) | 1 | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| E_h | Hartree energy, 2 hcR_inf | J | const.notation-atomic | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45–025010-46) |
| E_V | photon energy of the vacuum resonance (plasma and vacuum contributions to the dielectric tensor cancel) at given B and rho | keV | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| E_V | photon energy of the vacuum resonance (plasma and vacuum contributions to the dielectric tensor cancel) at given B and rho | keV | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| E_{ad} | adiabatic energy; resonant conversion occurs when E >~ E_ad | keV | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| E_{ce} | electron cyclotron energy hbar omega_ce = hbar e B/(m_e c) = 1158 B_14 keV | keV | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| E_{ci} | ion cyclotron energy hbar omega_ci = hbar Z e B/(A m_u c) = 0.635 (Z/A) B_14 keV | keV | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| E_{ci} | ion cyclotron energy hbar omega_ci = hbar Z e B/(A m_u c) = 0.635 (Z/A) B_14 keV | keV | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| F_\lambda | flux density per unit wavelength | erg cm^-2 s^-1 A^-1 | obs.notation-rxj1856-photometry | kaplan_2011_rxj1856_spectrum (Sec. 3.1 and Table 5; Sec. 2.3 and Table 2) |
| f_\nu | Eddington factor (second moment of P_ν over J_ν) | 1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| f_B | slowly varying function of B, f(B) ~ [3 delta_V/(q + m)]^(1/2) ~ 1 (Ho & Lai 2003, lines 444-457 of the local copy; VAL06 uses f_B ~ 1) | 1 | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| f_{D} | distance factor b_0/b of the epoch propagation | 1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| F_{ph} | local radial photon energy flux density at the stellar surface | erg cm^-2 s^-1 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| f_{T} | time factor of the epoch propagation with light-time effects | 1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| f_{V} | velocity factor of the epoch propagation with light-time effects | 1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| g | surface gravity (local proper) g = GM/(R^2 sqrt(1 - r_g/R)), taken uniform through the atmosphere | cm s^-2 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| g | surface gravity (local proper) g = GM/(R^2 sqrt(1 - r_g/R)), taken uniform through the atmosphere | cm s^-2 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| G'_{sRGB} | non-linear (encoded) sRGB green value | 1 | colour.notation-srgb-nonlinear | stokes_1996_srgb (eq. 1.2a, 1.2b) |
| G_{8bit} | 8-bit digital code value of the green channel | 1 (counts) | colour.notation-srgb-digital | stokes_1996_srgb (eq. 1.3-1.6) |
| G_{BP} | Gaia blue-photometer band magnitude | mag | obs.notation-gaia-photometry | riello_2021_gaia_photometry (Sect. 7, Eq. (7), p. 15 (refs/riello_2021_gaia_photometry.txt:931-940); Eqs. (13)-(17), p. 17; Table 3, p. 18) |
| G_{Gaia} | Gaia G-band magnitude | mag | obs.notation-gaia-photometry | riello_2021_gaia_photometry (Sect. 7, Eq. (7), p. 15 (refs/riello_2021_gaia_photometry.txt:931-940); Eqs. (13)-(17), p. 17; Table 3, p. 18) |
| g_{rad} | radiative acceleration | cm s^-2 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| G_{RP} | Gaia red-photometer band magnitude | mag | obs.notation-gaia-photometry | riello_2021_gaia_photometry (Sect. 7, Eq. (7), p. 15 (refs/riello_2021_gaia_photometry.txt:931-940); Eqs. (13)-(17), p. 17; Table 3, p. 18) |
| G_{sRGB} | linear sRGB green tristimulus value | 1 | colour.notation-srgb-linear | stokes_1996_srgb (eq. 1.1) |
| H_0 | σT_eff^4/(4π), the required Eddington flux integrated over frequency | erg s^-1 cm^-2 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| h_\nu | surface Eddington factor (first moment of P_ν over J_ν) | 1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| H_\nu | Eddington flux at column density m | erg s^-1 cm^-2 Hz^-1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| H_\rho | density scale height along the ray, \|ds/d ln ρ\| | cm | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| I^\infty | observed intensity at the observer (at infinity) corresponding to I_0 | erg s^-1 cm^-2 sr^-1 (inferred, as I_0) | gr.notation-radiation | beloborodov_2002_ray_tracing (Appendix, eq. 11 and text above it) |
| I_0 | surface radiation intensity at the emission point as a function of emission angle alpha (emitter frame); the paper relates it to the observed I through I/nu^4 = I_0/nu_0^4 | erg s^-1 cm^-2 sr^-1 (inferred frequency-integrated, from the I/nu^4 scaling) | gr.notation-radiation | beloborodov_2002_ray_tracing (Appendix, eq. 11 and text above it) |
| I_\nu(m,-\mu) | specific intensity of ingoing radiation (direction -μ) | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| I_\nu(m,\mu) | specific intensity of outgoing radiation at column depth y, frequency ν, direction μ | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| I_\nu^j | specific intensity in normal mode j (j = 1 X, 2 O) as a function of depth and direction; +mu outgoing, -mu ingoing | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| I_\nu^j | specific intensity in normal mode j (j = 1 X, 2 O) as a function of depth and direction; +mu outgoing, -mu ingoing | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| I_{1/2}(\chi,\tau) | standard Fermi-Dirac integral of order 1/2 for relativistic electrons (given by a fit in Chabrier & Potekhin 1998 per PY01) | 1 | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| j | normal-mode index, j = 1 extraordinary (X), j = 2 ordinary (O) | 1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| J_\nu | mean specific intensity, integral over 0..1 of P_nu dmu; the moment-equation J_nu equals the transfer-solution mean only as the temperature correction vanishes (McPHAC) | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| J_\nu^j | mean intensity of mode j, (1/2)∫_0^1 dμ [I(μ)+I(-μ)] | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| k^\mu | photon tangent four-vector dx^mu/dlambda_aff (normalisation source-specific: Bel02 u_t = 1; Pech83 eq. 2.2) | s cm^-1 for u^t; 1 for u^r; cm^-1 for u^psi (t in s, lambda in cm; source normalisation u_t = 1 is source-specific, conversion to the adopted metric convention unspecified) | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| k_\nu | total opacity k_ν = κ_ν + σ_T,r (called total absorption opacity in the paper) | cm^2 g^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| K_\parallel | thermal conductivity along the magnetic field B (GH83 eq. 6) | erg cm^-1 s^-1 K^-1 | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| K_\perp | thermal conductivity perpendicular to the magnetic field B (GH83 eq. 6); the anisotropic model assumes K_perp << K_par, the source also considers equality | erg cm^-1 s^-1 K^-1 | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| k_R | Rosseland mean total absorption opacity | cm^2 g^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| K_{ab}^{cond} | thermal conductivity tensor K_ab, spatial indices a, b in {x, z} in the horizontal/vertical frame (GH83 eqs. 7-8); not the mode indices i, j | erg cm^-1 s^-1 K^-1 | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| L_{ph}^\infty | apparent photon luminosity for a distant observer, L_ph^inf = (1 - r_g/R) L_ph | erg s^-1 | gr.notation-potekhin | potekhin_pons_page_2015_thermal (Sect. 2.1 (eqs. 1-3) and Sect. 5 (eqs. 39-40)) |
| L_{ph}^\infty | apparent photon luminosity for a distant observer, L_ph^inf = (1 - r_g/R) L_ph | erg s^-1 | surface.notation-observer-frame | potekhin_ho_chabrier_2016_review (Sect. 2.1, p. 4; Sect. 2.3, p. 8) |
| L_{slab} | depth of the isothermal reservoir below the surface (written L in GH83 Fig. 1) | cm | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| M | total gravitational stellar mass, M = M_R | g | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| M | total gravitational stellar mass, M = M_R | g | gr.notation-potekhin | potekhin_pons_page_2015_thermal (Sect. 2.1 (eqs. 1-3) and Sect. 5 (eqs. 39-40)) |
| m | column density (mass per unit area above the point), dm = -rho dz | g cm^-2 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| m | column density (mass per unit area above the point), dm = -rho dz | g cm^-2 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| M | total gravitational stellar mass, M = M_R | g | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| m_e | electron mass | kg | const.notation-particles | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45, 025010-47, 025010-49 (footnote f)) |
| m_n | neutron mass | kg | const.notation-particles | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45, 025010-47, 025010-49 (footnote f)) |
| m_p | proton mass | kg | const.notation-particles | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45, 025010-47, 025010-49 (footnote f)) |
| m_p/m_e | proton-electron mass ratio | 1 | const.notation-particles | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45, 025010-47, 025010-49 (footnote f)) |
| M_r | gravitational mass inside a sphere of radius r, M_r = 4 pi int_0^r rho r'^2 dr' | g | gr.notation-potekhin | potekhin_pons_page_2015_thermal (Sect. 2.1 (eqs. 1-3) and Sect. 5 (eqs. 39-40)) |
| m_u | atomic mass constant m_u (= unified atomic mass unit u) | kg | const.notation-particles | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45, 025010-47, 025010-49 (footnote f)) |
| M_{mod} | degree of modulation: (peak intensity - minimum intensity)/minimum intensity (GH83 eq. 23) | 1 | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| n_e | electron number density | cm^-3 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| N_H | interstellar hydrogen column density along the line of sight to the source | cm^-2 | obs.notation-rxj1856-spectral | ho_2007_rxj1856_thin_atm (Sec. 1; Sec. 2; Tables 1-2; Sec. 5.1) |
| N_H | interstellar hydrogen column density along the line of sight to the source | cm^-2 | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| n_L | Landau level index n = 0..n_max in the electron-pressure sum | 1 | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| n_{p,Gaia} | mean source flux in a Gaia band in electrons per second (e- s^-1), not normalised by telescope area; in the XP chapter n_p(lambda) is the photon flux distribution | s^-1 | obs.notation-gaia-photometry | riello_2021_gaia_photometry (Sect. 7, Eq. (7), p. 15 (refs/riello_2021_gaia_photometry.txt:931-940); Eqs. (13)-(17), p. 17; Table 3, p. 18) |
| N_{pix} | total number of pixels in the image (paper symbol N) | 1 | colour.notation-reinhard-parameters | reinhard_2002_tonemapping (p. 4, eq. 1-2 text) |
| P_\nu(m,\mu) | Feautrier variable: half the sum of outgoing and ingoing specific intensity, (I(y,μ)+I(y,-μ))/2 | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| P_g | gas pressure (electrons plus ions) | dyn cm^-2 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| P_g | gas pressure (electrons plus ions) | dyn cm^-2 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| P_g | gas pressure (electrons plus ions) | dyn cm^-2 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| P_r | relativistic unit of pressure m_e c^2/lambda_e^3 = 1.4218e25 dyn cm^-2 | dyn cm^-2 | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| P_{jump} | non-adiabatic jump probability; 1-P_jump is the O<->X conversion probability | 1 | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| P_{pupil} | telescope pupil area in the Gaia flux integral | cm^2 | obs.notation-gaia-photometry | riello_2021_gaia_photometry (Sect. 7, Eq. (7), p. 15 (refs/riello_2021_gaia_photometry.txt:931-940); Eqs. (13)-(17), p. 17; Table 3, p. 18) |
| P_{spin} | spin (rotation) period of the neutron star | s | obs.notation-rxj1856-timing | vankerkwijk_kaplan_2008_rxj1856_timing (Table 2 note) |
| P_{spin} | spin (rotation) period of the neutron star | s | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| q_z | upward heat-flux component at the surface, q_z(surface) = sigma T^4 (GH83 eqs. 10-11) | erg cm^-2 s^-1 | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| Q_{s} | cohesive energy of the condensed phase | eV | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| R | stellar (circumferential) radius | cm | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| r | radial coordinate measured from the stellar centre (areal radius) | cm | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| R | stellar (circumferential) radius | cm | gr.notation-potekhin | potekhin_pons_page_2015_thermal (Sect. 2.1 (eqs. 1-3) and Sect. 5 (eqs. 39-40)) |
| r | radial coordinate measured from the stellar centre (areal radius) | cm | gr.notation-potekhin | potekhin_pons_page_2015_thermal (Sect. 2.1 (eqs. 1-3) and Sect. 5 (eqs. 39-40)) |
| R | stellar (circumferential) radius | cm | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| R'_{sRGB} | non-linear (encoded) sRGB red value | 1 | colour.notation-srgb-nonlinear | stokes_1996_srgb (eq. 1.2a, 1.2b) |
| R^\infty | apparent radius of the star for a distant observer, R^inf = R/sqrt(1 - r_g/R) = R(1 + z_g) | cm | gr.notation-potekhin | potekhin_pons_page_2015_thermal (Sect. 2.1 (eqs. 1-3) and Sect. 5 (eqs. 39-40)) |
| R^\infty | apparent radius of the star for a distant observer, R^inf = R/sqrt(1 - r_g/R) = R(1 + z_g) | cm | surface.notation-observer-frame | potekhin_ho_chabrier_2016_review (Sect. 2.1, p. 4; Sect. 2.3, p. 8) |
| R^\infty_{bb}/D | X-ray blackbody normalisation (apparent radius per distance) | km kpc^-1 | obs.notation-rxj1856-photometry | kaplan_2011_rxj1856_spectrum (Sec. 3.1 and Table 5; Sec. 2.3 and Table 2) |
| r_e | classical electron radius, alpha^2 a_0 | m | const.notation-atomic | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45–025010-46) |
| r_g | Schwarzschild radius 2GM/c^2 of the star | cm | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| r_g | Schwarzschild radius 2GM/c^2 of the star | cm | gr.notation-potekhin | potekhin_pons_page_2015_thermal (Sect. 2.1 (eqs. 1-3) and Sect. 5 (eqs. 39-40)) |
| r_g | Schwarzschild radius 2GM/c^2 of the star | cm | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| R_{8bit} | 8-bit digital code value of the red channel | 1 (counts) | colour.notation-srgb-digital | stokes_1996_srgb (eq. 1.3-1.6) |
| R_{cap} | radius of the polar cap, R_cap = R^(3/2) [2 pi/(c P_spin)]^(1/2) << R (GH83 eq. 1, where the period is written p) | cm | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| R_{em} | physical size of the emission region | km | obs.notation-rxj1856-spectral | ho_2007_rxj1856_thin_atm (Sec. 1; Sec. 2; Tables 1-2; Sec. 5.1) |
| R_{em}^\infty | apparent (radiation) radius at infinity, R^inf = Rem (1+z_g) | km | obs.notation-rxj1856-spectral | ho_2007_rxj1856_thin_atm (Sec. 1; Sec. 2; Tables 1-2; Sec. 5.1) |
| R_{Ry} | Rydberg constant | m^-1 | const.notation-atomic | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45–025010-46) |
| R_{sRGB} | linear sRGB red tristimulus value | 1 | colour.notation-srgb-linear | stokes_1996_srgb (eq. 1.1) |
| S_\nu(m,\mu) | source function | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| S_\nu^j | source function of normal mode j | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| S_\nu^j | source function of normal mode j | erg s^-1 cm^-2 Hz^-1 sr^-1 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| S_{Gaia} | Gaia system overall response function, including the scaling from photons to electrons (passband) | e- photon^-1 | obs.notation-gaia-photometry | riello_2021_gaia_photometry (Sect. 7, Eq. (7), p. 15 (refs/riello_2021_gaia_photometry.txt:931-940); Eqs. (13)-(17), p. 17; Table 3, p. 18) |
| S_{stim} | stimulated-emission factor S_stim = 1 - exp(-E/k_B T) | 1 | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| S_{vis} | visible area of the stellar surface | cm^2 | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| t | Schwarzschild coordinate time | s | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| T^{(0)} | present best-guess temperature profile T_0(y) of the iteration, T = T_0 + Delta T | K | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| T_* | temperature from the Appendix-B Ts fits evaluated without the neutrino-limit correction (B.4), used to define the flux F_b = sigma_SB T_*^4 (PPP15 App. B item 4) | K | surface.notation-temperatures | potekhin_pons_page_2015_thermal (eqs. 39, 40; Sect. 5.1 (L1942-1943); App. B (L2674)) |
| T_6 | T / (10^6 K) | 1 | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| T_b | temperature at the bottom of the heat-blanketing envelope (crust/envelope boundary, density rho_b) (PPP15 Sect. 5.1) | K | surface.notation-temperatures | potekhin_pons_page_2015_thermal (eqs. 39, 40; Sect. 5.1 (L1942-1943); App. B (L2674)) |
| T_r | relativistic unit of temperature m_e c^2/k_B = 5.930e9 K | K | eos.notation-electron-gas | potekhin_yakovlev_2001_cooling (p. 1 (b, lambda_e, omega_c) and p. 4, eq. 12 and text) |
| T_s | local effective surface temperature in the local reference frame at the stellar surface | K | rt.notation-mcphac | haakonsen_2012_mcphac (Sects. 2–3.4, eqs. 1–7, 10, 19–22; Appendix A, eq. A4) |
| T_s | local effective surface temperature in the local reference frame at the stellar surface | K | surface.notation-temperatures | potekhin_pons_page_2015_thermal (eqs. 39, 40; Sect. 5.1 (L1942-1943); App. B (L2674)) |
| T_{bg} | background temperature of the stellar surface outside the cap (written T in GH83 eqs. 2-3) | K | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| T_{cap} | temperature of the heated polar cap (GH83 eq. 2-3) | K | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| T_{crit} | critical temperature below which condensation occurs | K | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| T_{eff} | overall effective temperature of the star in the local frame, 4 pi sigma_SB R^2 T_eff^4 = L_ph = integral of F_ph dSigma (PPP15 eq. 39) | K | surface.notation-temperatures | potekhin_pons_page_2015_thermal (eqs. 39, 40; Sect. 5.1 (L1942-1943); App. B (L2674)) |
| T_{eff}^\infty | redshifted (apparent) effective temperature seen by a distant observer, T_eff^inf = T_eff sqrt(1 - r_g/R) = T_eff/(1 + z_g) | K | gr.notation-potekhin | potekhin_pons_page_2015_thermal (Sect. 2.1 (eqs. 1-3) and Sect. 5 (eqs. 39-40)) |
| T_{eff}^\infty | redshifted (apparent) effective temperature seen by a distant observer, T_eff^inf = T_eff sqrt(1 - r_g/R) = T_eff/(1 + z_g) | K | obs.notation-rxj1856-spectral | ho_2007_rxj1856_thin_atm (Sec. 1; Sec. 2; Tables 1-2; Sec. 5.1) |
| T_{eff}^\infty | redshifted (apparent) effective temperature seen by a distant observer, T_eff^inf = T_eff sqrt(1 - r_g/R) = T_eff/(1 + z_g) | K | surface.notation-temperatures | potekhin_pons_page_2015_thermal (eqs. 39, 40; Sect. 5.1 (L1942-1943); App. B (L2674)) |
| t_{ep} | reference epoch chosen by convention | s | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| T_{eq} | surface temperature at the magnetic equator (PPP15 App. B, T_eq) | K | surface.notation-temperatures | potekhin_pons_page_2015_thermal (eqs. 39, 40; Sect. 5.1 (L1942-1943); App. B (L2674)) |
| T_{fit} | best-fit blackbody temperature of a synthesised spectrum (GH83 Sect. III c; the source describes it in prose, T_fit is project notation) | K | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| T_{pole} | surface temperature at the magnetic pole (T_pole in GH83 eq. 14; T_p in PPP15 App. B) | K | surface.notation-temperatures | potekhin_pons_page_2015_thermal (eqs. 39, 40; Sect. 5.1 (L1942-1943); App. B (L2674)) |
| t_{pr} | time elapsed since the reference epoch over which the astrometric parameters are propagated (Julian years of 365.25 x 86400 s, TCB) | s | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| T_{res} | temperature of the isothermal heat reservoir below the surface slab (written T_0 in GH83 Fig. 1) | K | surface.notation-conduction | greenstein_hartke_1983_surface_T (eqs. 1-3, 6-11, 23; Sect. III c, pp. 289-290) |
| u_\nu^j | specific energy density of mode j, (2π/c)∫_-1^1 dμ′ I_ν^j(μ′) | erg cm^-3 Hz^-1 | rt.notation-val | van_adelsberg_lai_2006_magnetic_atm (Sects. 2-3.1, eqs. 6-9, 30-34) |
| u_r | relative standard uncertainty | 1 | const.notation-atomic | tiesinga_2021_codata2018 (Table XXXI, pp. 025010-45–025010-46) |
| u_{e} | (E_Be/E)^2, squared ratio of electron cyclotron energy to photon energy | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| u_{i} | (E_Bi/E)^2, squared ratio of ion cyclotron energy to photon energy | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| v_r | radial velocity of a source relative to the solar system barycentre (sixth kinematic parameter of the standard model) | cm s^-1 | obs.notation-gaia-sky | lindegren_2021_gaia_astrometry (Sect. 3.1, p. 6 (refs/lindegren_2021_gaia_astrometry.txt:432-437,470-475)) |
| v_{e} | (E_pe/E)^2, squared ratio of electron plasma energy to photon energy | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| v_{t} | apparent tangential velocity of the source, v_t = mu b (BL Sect. 4, apparent quantities); the true tangential velocity is a different quantity | cm s^-1 | obs.notation-gaia-propagation | butkevich_lindegren_2014_propagation (Sect. 4 (astrometric parameters), Sect. 5, Appendix D (tex l.517)) |
| X | CIE 1931 tristimulus value X | same as Y | colour.notation-cie-xyz | stokes_1996_srgb (eq. 1.1, 1.2) |
| x | CIE chromaticity coordinate x | 1 | colour.notation-chromaticity | stokes_1996_srgb (Table 0.2) |
| x_g | compactness r_g/R = 2GM/(c^2 R) at the stellar surface | 1 | gr.notation-geometry | beloborodov_2002_ray_tracing (Sect. 1-3 and Appendix) |
| x_g | compactness r_g/R = 2GM/(c^2 R) at the stellar surface | 1 | surface.notation-observer-frame | potekhin_ho_chabrier_2016_review (Sect. 2.1, p. 4; Sect. 2.3, p. 8) |
| x_{pix} | pixel x coordinate in the image | 1 | colour.notation-reinhard-parameters | reinhard_2002_tonemapping (p. 4, eq. 1-2 text) |
| Y | CIE 1931 tristimulus value Y (luminance-like; photometric normalisation not stated in this source) | source-unspecified | colour.notation-cie-xyz | stokes_1996_srgb (eq. 1.1, 1.2) |
| y | CIE chromaticity coordinate y | 1 | colour.notation-chromaticity | stokes_1996_srgb (Table 0.2) |
| Y_e | electron fraction Z/A | 1 | rt.notation-vacuum | van_adelsberg_lai_2006_magnetic_atm (Sect. 2, eqs. 1-5) |
| y_H | hydrogen column depth of the thin atmosphere (mass per area) | g cm^-2 | obs.notation-rxj1856-spectral | ho_2007_rxj1856_thin_atm (Sec. 1; Sec. 2; Tables 1-2; Sec. 5.1) |
| y_{pix} | pixel y coordinate in the image | 1 | colour.notation-reinhard-parameters | reinhard_2002_tonemapping (p. 4, eq. 1-2 text) |
| Z | CIE 1931 tristimulus value Z | same as Y | colour.notation-cie-xyz | stokes_1996_srgb (eq. 1.1, 1.2) |
| z | CIE chromaticity coordinate z | 1 | colour.notation-chromaticity | stokes_1996_srgb (Table 0.2) |
| z_g | gravitational redshift, 1 + z_g = (1 - r_g/R)^(-1/2) = (1 - x_g)^(-1/2) | 1 | obs.notation-rxj1856-spectral | ho_2007_rxj1856_thin_atm (Sec. 1; Sec. 2; Tables 1-2; Sec. 5.1) |
| z_g | gravitational redshift, 1 + z_g = (1 - r_g/R)^(-1/2) = (1 - x_g)^(-1/2) | 1 | surface.notation-observer-frame | potekhin_ho_chabrier_2016_review (Sect. 2.1, p. 4; Sect. 2.3, p. 8) |
| Z_{ion} | ion charge number Z (nuclear charge of a fully ionised ion) | 1 | rt.notation-cond-matter-dielectric | van_adelsberg_2005_condensed_surface (Sect. 2 (p. 2-4)) |
| z_{vert} | vertical distance z in the atmosphere (SPW09) | cm | rt.notation-spw09 | suleimanov_potekhin_werner_2009_mag_atm (Sect. 2, eqs. 1-25) |
| ZP_{AB} | Gaia passband zero point in the AB system, per band | mag | obs.notation-gaia-photometry | riello_2021_gaia_photometry (Sect. 7, Eq. (7), p. 15 (refs/riello_2021_gaia_photometry.txt:931-940); Eqs. (13)-(17), p. 17; Table 3, p. 18) |
| ZP_{VEG} | Gaia passband zero point in the VEGAMAG system, per band | mag | obs.notation-gaia-photometry | riello_2021_gaia_photometry (Sect. 7, Eq. (7), p. 15 (refs/riello_2021_gaia_photometry.txt:931-940); Eqs. (13)-(17), p. 17; Table 3, p. 18) |
