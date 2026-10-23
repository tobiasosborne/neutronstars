# Notation additions: camera lane (GAPS 1.4, 1.6)

Proposed rows for `truth/NOTATION_DECISIONS.md`, section "GR exterior, ray geometry and stellar
structure" (or a new "Camera frame" subsection). Same columns as the dictionary. No notation
entries exist for these yet, so "adopted" names the equation entry in `truth/entries/camera.toml`
that introduces the symbol. These rows were chosen to avoid collisions: bare `\beta` is the bending
angle, bare `g` is surface gravity, bare `n` and `u` are retired, `D` is distance, `\Lambda` is the
metric function, and `N_A`, `N_H` and `N_{pix}` are taken.

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `e^{st}_{(t)}` | time leg of the static orthonormal Schwarzschild tetrad, e = [c (1 - r_g/r)^(1/2)]^(-1) d/dt (u_obs/c) | s cm^-1 (component along d/dt) | `e_(t)` (gr.static-tetrad, MG09 eq. 2.2.6); FIDO time leg (James+15, a = 0) | gr.static-tetrad | |
| `e^{st}_{(r)}` | radial leg of the static tetrad, (1 - r_g/r)^(1/2) d/dr | 1 | `e_(r)` (gr.static-tetrad); `e_r-hat` (James+15 eq. A.3, gr.camera-aberration); `e_lambda-hat` (WS83, gr.static-tetrad-field-projection) | gr.static-tetrad | |
| `e^{st}_{(\theta)}` | polar leg of the static tetrad, r^(-1) d/dtheta_Sch | cm^-1 | `e_(vartheta)` (gr.static-tetrad); `e_theta-hat` (James+15 eq. A.3) | gr.static-tetrad | Label refers to `\theta_{Sch}`, not the RT ray angle `\theta`. |
| `e^{st}_{(\phi)}` | azimuthal leg of the static tetrad, (r sin theta_Sch)^(-1) d/dphi | cm^-1 | `e_(varphi)` (gr.static-tetrad); `e_phi-hat` (James+15 eq. A.3) | gr.static-tetrad | Schwarzschild azimuth, not the surface-element azimuth `\phi`. |
| `e_{(a)}` | camera orthonormal tetrad leg, a in {0, 1, 2, 3}; e_(0) = u_cam/c | e_(0)^t in s cm^-1; spatial legs as the static legs they reduce to | `e_(i)`, `e-hat_(i)` (gr.tetrad-lorentz-transformation, MG09 eq. 1.4.6); `e_x, e_y, e_z` (James+15 Sect. 2.1, gr.camera-aberration) | gr.tetrad-lorentz-transformation | Already the pipeline spelling (dictionary "Pipeline-level spellings"); this row makes it canonical. |
| `e^{(a)}` | dual (co)tetrad, e^(a)_mu e_(b)^mu = delta^a_b | inverse of e_(a) | `theta^(i)` (MG09 eqs. 1.4.1, 1.4.19, gr.tetrad-vector-components) | gr.tetrad-vector-components | MG09's theta^(i) is renamed: `\theta` is the RT ray angle. |
| `u^\mu_{cam}` | four-velocity of the camera, u = c gamma_cam (e_(0)^st + beta_cam v-hat_cam) | as `u^\mu_{obs}` (u^t dimensionless with t in s; u^r cm s^-1; u^theta, u^phi s^-1) | `u` (MG09 eq. 1.4.18, gr.camera-four-velocity-static-frame) | gr.camera-four-velocity-static-frame | Listed in the pipeline-spelling table without a row. |
| `\beta_{cam}` | camera speed relative to the local static observer, in units of c | 1 | `beta` (MG09 eq. 1.4.7; James+15 eq. A.9, speed relative to the FIDO) | gr.camera-four-velocity-static-frame | Not the bending angle `\beta`. |
| `\gamma_{cam}` | Lorentz factor of the camera relative to the local static observer, (1 - beta_cam^2)^(-1/2) | 1 | `gamma` (MG09 eq. 1.4.7) | gr.tetrad-lorentz-transformation | Not `\gamma_B` or `\gamma^2`. |
| `\hat v_{cam}` | unit 3-vector of the camera's direction of motion in the local static frame | 1 | `n^a`, `n` (MG09 eqs. 1.4.7, 1.4.18); `B` (James+15 Sect. 2.1, e_y = B) | gr.tetrad-lorentz-transformation | MG09 n is renamed (bare n is retired); James's B is not a magnetic field. |
| `\Lambda_{boost}` | Lorentz boost matrix between the static and camera tetrads | 1 | `Lambda^a_b`, `A_i^k` (MG09 eqs. 1.4.6-1.4.7) | gr.tetrad-lorentz-transformation | Not the metric function `\Lambda`. Sign convention: see gr.tetrad-lorentz-transformation notes. |
| `\hat N_{cam}` | unit vector in the camera frame pointing from the camera toward the apparent source of an incoming ray (opposite to propagation) | 1 | `N`, `N_x, N_y, N_z` (James+15 eq. A.8, gr.camera-aberration) | gr.camera-aberration | Photon propagation direction in the camera frame is -N_cam. |
| `\hat n_{st}` | propagation direction of an incoming ray measured by the local static observer | 1 | `n_F` (James+15 eq. A.9) | gr.camera-aberration | |
| `\theta_{cam}` | polar angle on the camera's local sky, measured from camera axis e_(3) (James e_z) | rad | `theta_cs` (James+15 eq. A.8) | gr.camera-aberration | |
| `\phi_{cam}` | azimuth on the camera's local sky about e_(3), from e_(1) (James e_x) | rad | `phi_cs` (James+15 eq. A.8) | gr.camera-aberration | |
| `\nu_{cam}` | photon frequency measured in the camera frame | Hz | `f_c` (James+15 eq. A.16); `nu_0` (DA09 eq. 70, observed) | gr.frequency-ratio-general | Listed in the pipeline-spelling table without a row. |
| `\nu_{st}` | photon frequency measured by the local static observer at the camera's position | Hz | text/other: FIDO frequency (James+15 eq. A.11, E_F) | gr.doppler-factor-sr | |
| `g_{shift}` | frequency ratio nu_cam/nu = (p.u_cam)/(p.u_em), with nu the emitter-frame frequency | 1 | `g` (DA09 eq. 70, g = nu_0/nu); `f_c/f'` (James+15 eq. A.16, emitter at rest at infinity) | gr.frequency-ratio-general | Already the pipeline spelling; bare g is surface gravity. Static emitter at R, camera at infinity: g_shift = (1 - x_g)^(1/2) = 1/(1 + z_g). |
| `\mathcal{D}_{cam}` | special-relativistic Doppler factor nu_cam/nu_st = (1 - beta_cam^2)^(1/2)/(1 - beta_cam N_y) | 1 | first factor of `f_c/f'` (James+15 eq. A.16) | gr.doppler-factor-sr | `D` is distance. |
| `f_{cam}` | camera focal length in pixel units (projection constant) | 1 (pixels) | none (project definition) | gr.camera-projection-convention | Requires f_cam > 0. Not `f_\nu`, `f_{col}`, `f_B`. |
| `r_{pix}` | radial distance of a pixel from the optical-axis pixel | 1 (pixels) | none (project definition) | gr.camera-projection-convention | |
| `\theta_{field}` | field angle between N_cam and the optical axis a_cam | rad | none (project definition) | gr.camera-projection-convention | Pinhole: tan theta_field = r_pix/f_cam; equidistant: theta_field = r_pix/f_cam (open PRD decision). Requires f_cam > 0 and 0 <= theta_field <= pi; the full field of view, not theta_field, may exceed pi. |
| `\phi_{field}` | image-plane azimuth of a pixel about the optical axis, from image-right | rad | none (project definition) | gr.camera-projection-convention | |
| `\hat a_{cam}` | optical-axis (pointing) unit vector in the camera frame | 1 | none (project definition; PRD The shot, Pointing) | gr.camera-projection-convention | |

Row count: 24.

## Bridges (camera review r1 cross-entry remarks)

| pair | distinction |
|---|---|
| lapse `\alpha` vs emission angle `\alpha` | James+15's lapse `alpha` (eq. A.16, gr.camera-blueshift-celestial) is a metric function evaluated at the camera, not an emission angle; keep it spelled `alpha` (lapse) only inside quoted James equations. |
| axial `b` vs total impact parameter `b` | James+15's `b` is the ray's axial angular momentum p_phi with p_t = -1; the RT total impact parameter is a different quantity (they coincide only for equatorial rays). |
| direction amplitude `\psi` vs orbital angle `\psi` | MG09's `psi` (upsilon = upsilon^(0) e_(0) + psi n) is the spatial-direction amplitude (psi = 1 for null); it is not an orbital angle. |
| direction angle `\chi` vs plasma susceptibility `\chi` | MG09's polar direction angle `chi` (gr.tetrad-vector-components, gr.camera-four-velocity-static-frame) is not a plasma susceptibility. |
| physical `p^\mu` vs normalized tangent | Physical photon momentum p^mu carries energy h nu; the normalized ray tangent has p_t = -1 (James+15). Related by p^(0) = E_obs/c = -p_mu u^mu/c. |
| observed `\nu_0` vs emitter frequency | Dexter & Agol's `nu_0` is the observed frequency (g = nu_0/nu, nu emitted); Beloborodov's spelling assigns the subscript to the emitter frequency, i.e. the reverse. Project spelling: `\nu_{cam}` observed, `\nu` emitted. |
