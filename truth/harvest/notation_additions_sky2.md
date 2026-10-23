# Notation additions, lane sky2 (all from obs.notation-gaia-propagation)

Dictionary row format: symbol | meaning | units | source spellings (entries) | adopted | note. Reuses (does not redefine) the rows of notation_additions_sky.md: `\alpha_{RA}`, `\delta_{Dec}`, `\varpi_{plx}`, `\mu_{\alpha*}`, `\mu_\delta`, `v_r`, `A_{u}`.

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `t_{ep}` | reference epoch chosen by convention | s | `T_0` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | sky2-r1 fix. Mission-midpoint preference (Lindegren 2012) and J2016.0 (Lindegren 2021) are separate claims, not part of the meaning. |
| `t_{pr}` | time elapsed since the reference epoch over which the astrometric parameters are propagated (Julian years of 365.25 x 86400 s, TCB) | s (catalogue: yr) | `t = T - T_0` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\vec{u}_{bary}` | barycentric coordinate direction to the source, unit vector b/|b| in the ICRS | 1 | `u` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\vec{r}_{tri}` | fixed unit vector coinciding with the barycentric direction at the reference epoch; third axis of the normal triad | 1 | `r` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\vec{p}_{tri}` | unit vector in the direction of increasing right ascension (first axis of the normal triad) | 1 | `p` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\vec{q}_{tri}` | unit vector in the direction of increasing declination (second axis of the normal triad) | 1 | `q` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\vec{\mu}_{pm}` | proper motion vector du/dt of the barycentric direction | rad s^-1 (catalogue: mas yr^-1) | `mu (vector)` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\mu_{r}` | radial proper motion mu_r = v_r varpi / A (radial velocity in the unit of proper motion), sixth astrometric parameter | rad s^-1 (catalogue: mas yr^-1) | `mu_r` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\mu_{tot}` | total tangential proper motion, (mu_alpha*^2 + mu_delta^2)^(1/2) = v_t / b | rad s^-1 (catalogue: mas yr^-1) | `mu` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\vec{b}_{bary}` | barycentric position vector of the source: the apparent position b_app(T_obs) = b_true(T_em) from BL Sect. 4 onward, where subscripts are omitted; the true position is a different quantity | cm | `b (vector)` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `b_{bary}` | barycentric distance of the source, varpi = A/b | cm | `b = |b|` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\vec{v}_{bary}` | space velocity of the source relative to the solar system barycentre: constant TRUE velocity in the uniform rectilinear model (BL Sect. 3, eq. 1); APPARENT velocity from BL Sect. 4 onward; the two differ by the light-time factor 1/(1 - v_r^app/c) of BL eq. (19) | cm s^-1 | `v (vector)` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `v_{t}` | apparent tangential velocity of the source, v_t = mu b (BL Sect. 4); the true tangential velocity is a different quantity | cm s^-1 | `v_t` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `A_{vel}` | astronomical unit expressed in km yr s^-1: A_V = 4.740 470 446 | km yr s^-1 | `A_V` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\tau_{A}` | light-travel time of one astronomical unit, tau_A = A/c = 499.004 784 s | s | `tau_A` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\tau_{lt}` | light-travel time T_obs - T_em from the source to the barycentre; tau_0 = b_0/c at the reference epoch | s | `tau` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `f_{D}` | distance factor b_0/b of the epoch propagation | 1 | `f_D` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `f_{T}` | time factor of the epoch propagation with light-time effects | 1 | `f_T` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `f_{V}` | velocity factor of the epoch propagation with light-time effects | 1 | `f_V` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\varepsilon_{T}` | small parameter f_T - 1 of the light-time expansion | 1 | `epsilon_T` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\varepsilon_{V}` | small parameter f_V - 1 of the light-time expansion | 1 | `epsilon_V` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
| `\vec{v}_{true}` | true (constant) space velocity of the uniform rectilinear model, BL eq. (1); distinct from `\vec{v}_{app}` | cm s^-1 | `v (true)` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | sky2-r1 cross-entry fix. |
| `\vec{v}_{app}` | apparent space velocity through the astrometric parameters, v = (A_V/varpi)(r mu_r + mu), BL eq. (30); distinct from `\vec{v}_{true}` | cm s^-1 | `v (apparent)` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | sky2-r1 cross-entry fix. |
| `\mu_{tot,0}` | total tangential proper motion at the reference epoch, mu_0 = (mu_alpha*0^2 + mu_delta0^2)^(1/2) of BL | rad s^-1 (catalogue: mas yr^-1) | `mu_0` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | COLLISION: bare `\mu_0` is the vacuum magnetic permeability (const.notation-fundamental, `μ_0`; see NOTATION_DECISIONS.md) and `\mu` the ray cosine. BL's mu_0 is mapped to `\mu_{tot,0}`; the sky2 latex fields keep BL's spelling verbatim and must be read through this mapping. |
| `\Delta\theta_{lt}` | angular-position error caused by neglecting light-time effects in the epoch propagation (first-order estimate) | rad | `Delta theta` (obs.notation-gaia-propagation) | obs.notation-gaia-propagation | |
