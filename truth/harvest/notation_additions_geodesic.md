# Notation additions, geodesic lane (GAPS 1.2, 1.3, 1.7)

Proposed rows for `truth/NOTATION_DECISIONS.md`. Same columns as the Dictionary, plus the
section each row belongs in. Entries cited are in `truth/entries/gr_geodesic.toml`
(all `proposed`). The entry `latex` fields keep the source spellings (Dexter & Agol 2009 =
DA09, Beloborodov 2002 = Bel02, Pechenick 1983 = Pech83). These rows are the bridge to
the project spellings.

## Choice of spelling for the conserved specific energy and angular momentum

The lead suggested `E_{geo}`, `L_{geo}`. This file proposes `\tilde E` and `\tilde L` instead,
for two reasons:

1. The subscript `geo` already means "geometric units" in the dictionary (`u_{geo}` =
   GM/(c^2 r), `u_{b,geo}` = GM/(c^2 b)). So `E_{geo}` would read as "E in geometric units".
   That is wrong in CGS, where the quantity is the dimensionless (1 - r_g/r) dt/dtau, or
   (cm^2 s^-1) for L.
2. The canonical "Pipeline-level spellings" table already reserves `\tilde E_{cam}` and
   `\tilde L_{cam}` for the camera's specific energy and angular momentum. With
   `\tilde E`, `\tilde L` as the generic geodesic constants, the camera values are just
   these symbols evaluated on the camera worldline.

Neither spelling collides literally with an existing `sym`. Bare `E` is photon energy, and
`\tilde E_{C}`, `\tilde E_{pe}` are condensed-surface energies; these are distinct strings.
Bare `L` stays retired. If the lead prefers `E_{geo}`, `L_{geo}`, substitute them
one-for-one. Nothing else in this file depends on the choice.

## Proposed rows

Section "GR exterior, ray geometry and stellar structure":

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `\tilde E` | specific energy (per unit rest mass) of a timelike geodesic, conserved: -u_t in geometric units, (1 - r_g/r) dt/dtau in CGS | 1 | `E` (DA09 eqs. 6-10; gr.timelike-constants-schwarzschild, gr.timelike-radial-equation-schwarzschild, gr.timelike-dt-dtau-schwarzschild, gr.timelike-periastron-cubic, gr.timelike-angular-momentum-from-periastron, gr.timelike-unbound-condition, gr.timelike-speed-at-infinity, gr.timelike-proper-time-quadrature) | gr.timelike-constants-schwarzschild | Per unit rest mass: follows from the derived norm -delta_1 = -1 of dx/dlambda. DA09 Sect. 3 also writes gamma = E/m with explicit m; the two agree only for m = 1. Camera value: `\tilde E_{cam}`. Hyperbolic iff `\tilde E` > 1. |
| `\tilde L` | specific angular momentum (per unit rest mass) of a timelike geodesic about the orbit-plane normal, conserved: u_phi = r^2 dphi/dtau | cm^2 s^-1 (geometric: cm) | `L_z` (DA09; signed component; \|L_z\| is the total angular-momentum magnitude for the equatorial orbit, Carter Q = 0) | gr.timelike-constants-schwarzschild | Geometric L_z = `\tilde L`/c. Camera value: `\tilde L_{cam}`. Sign = orbit orientation, a free choice. L and L_z are signed synonyms only after the equatorial plane is chosen. |
| `\tau_{prop}` | proper time along a timelike worldline, dtau^2 = -ds^2/c^2 (signature -,+,+,+) | s | `\lambda` (DA09 affine parameter; for delta_1 = 1, lambda = tau + C after choosing orientation; geometric units, a length = c tau) | gr.timelike-constants-schwarzschild | Bare tau is taken (optical depths `\tau_\nu`, `\tau_T`). On the camera worldline this is `\tau_{cam}`. DA09 never names proper time; Mueller & Grave 2009 (refs/muller_grave_2009_spacetimes.tex line 700) print 'For timelike geodesics, however, we replace the affine parameter by the proper time $\tau$.' The explicit dtau^2 = -ds^2/c^2 is the timelike constraint there with kappa = -1 for signature +2 (tex line 707). |
| `r_{peri}` | periastron (minimum Schwarzschild radial coordinate) of an unbound timelike geodesic: the largest simple root of `\mathcal{R}`(r) with d`\mathcal{R}`/dr > 0 | cm | `r_p` (lane brief, PRD "Periastron distance") | gr.timelike-periastron-cubic | Not `r_p`: subscript p already means proton (`m_p`) and polar (`B_p`). |
| `r_{cam}` | Schwarzschild radial coordinate of the camera | cm | `r_c` (lane brief, pipeline.null-geodesic-trace) | gr.bending-integral-finite-radius | Spelling already reserved in "Pipeline-level spellings"; this row makes it a dictionary row. Distinct from `r_0` (Pech83 static observer, r_0 -> infinity) and from `D`. |
| `r_{tp}` | radial turning point (periastron) of a photon orbit with impact parameter b: largest root of r^3 - b^2 r + b^2 r_g = 0 | cm | `u_*` (DA09 Sect. 3, as u = 1/r with M = 1) | gr.bending-integral-finite-radius | Exists outside R only for b > b_max (R > 3GM/c^2). |
| `v^\infty` | asymptotic speed of an unbound timelike geodesic relative to static observers at infinity, c (1 - `\tilde E`^-2)^(1/2) | cm s^-1 | (none printed; derived) | gr.timelike-speed-at-infinity | Superscript infinity per the observer-frame convention. |
| `\phi_{Sch}` | Schwarzschild (Boyer-Lindquist at a = 0) azimuthal coordinate | rad | `\phi` (DA09 eqs. 1, 9) | gr.schwarzschild-line-element | Bare `\phi` is the surface-element azimuth about the line of sight (Bel02). For an equatorial orbit, `\phi_{Sch}` = `\phi_{orb}` up to the choice of origin. Pairs with `\theta_{Sch}`. |
| `\mathcal{R}(r)` | radial function of the separated geodesic equations (DA09 eq. 6); for a = 0, delta_1 = 1, Q = 0: r^4 E^2 - (r^2 - 2Mr)(L^2 + r^2) = r^4 (dr/dtau)^2 | cm^4 (geometric) | `R` (DA09 eq. 6) | gr.timelike-radial-equation-schwarzschild | Bare `R` is the stellar radius. |
| `\psi_{cam}` | finite swept polar angle of a photon between the stellar surface and a camera at r_cam (camera-referenced), int K_b dr with K_b(r) = b/[r^2 sqrt(1 - b^2 (1 - r_g/r)/r^2)] | rad | `psi` (Bel02 eq. 9 with upper limit r_cam; Pech83 eq. 2.12 with lower limit M/r_cam) | gr.bending-integral-finite-radius | Distinct from the dictionary `\psi` (infinity-referenced, psi = 0 at infinity), to which it tends as r_cam -> infinity. |
| `\delta_1` | geodesic type switch: 0 null, 1 timelike (DA09) | 1 | `\delta_1` (DA09 eqs. 6-7) | gr.timelike-constants-schwarzschild | Distinct from `\delta_{L}` (Reinhard offset). |

Section "GR exterior" (Kerr parent only; outside PRD scope; needed so that the a = 0 reduction
can be read against its printed parent):

| symbol | meaning | units | source spellings (entries) | adopted | note |
|---|---|---|---|---|---|
| `a_{spin}` | Kerr parameter; DA09: "the angular momentum of the black hole" in G = c = M = 1 | cm (geometric) | `a` (DA09 eqs. 1-3) | gr.kerr-line-element-bl | Bare a would collide with `a_{key}`, `a_q`, `a_0`, `a^{B.2}`. Set to 0 everywhere in this project. |
| `\Delta_K` | Kerr function r^2 - 2Mr + a^2 (DA09 eq. 2 with M = 1) | cm^2 | `\Delta` (DA09 eq. 2) | gr.kerr-line-element-bl | a = 0: r^2 - r_g r. Kept apart from `\Delta_\rho`, `\Delta T`. |
| `\rho_K^2` | r^2 + a^2 cos^2 theta (DA09 eq. 2) | cm^2 | `\rho^2` (DA09 eq. 2) | gr.kerr-line-element-bl | Bare `\rho` is mass density. a = 0: r^2. |
| `\Sigma_K^2` | (r^2 + a^2)^2 - a^2 Delta sin^2 theta (DA09 eq. 3) | cm^4 | `\Sigma^2` (DA09 eq. 3) | gr.kerr-line-element-bl | a = 0: r^4. |
| `\Theta_K(\theta)` | polar function of the separated geodesic equations (DA09 eq. 7) | cm^2 (geometric) | `\Theta` (DA09 eq. 7) | gr.kerr-line-element-bl | Avoids the pipeline's old `Theta` (field-normal angle, now `\theta_{Bn}`). |
| `\mathcal{Q}` | Carter constant (DA09) | cm^2 (geometric) | `\mathcal{Q}` (DA09 eqs. 6-7) | gr.timelike-constants-schwarzschild | 0 for the equatorial orbit chosen at a = 0. Distinct from `Q_\nu` (Stokes) and `Q_{s}` (cohesive energy). |

## Changed rows (aliases only)

| symbol | change |
|---|---|
| `\lambda_{aff}` | add source spelling `lambda` (DA09 eq. 4 Hamilton-Jacobi form, eqs. 8-10). Normalisation: for delta_1 = 1, lambda = tau + C after choosing orientation (derived in gr.timelike-constants-schwarzschild; geometric lambda_geom = c tau in CGS). For delta_1 = 0 it is fixed by the separated action, source-specific like the Bel02 and Pech83 ones. The r1 "disputed" note stands. |
| `GM/c^2` | add source spelling `M (G=c=M=1)` (DA09, which sets M = 1). The geometric-unit mass in the a = 0 reduction is this symbol. |
| `r_g` | add note: g_rr = (1 - r_g/r)^-1 in the exterior, now carried by gr.schwarzschild-line-element (replaces the rejected PPP15-based metric rows). |
| `\Lambda` | add note: exterior identity e^(-2 Lambda) = 1 - r_g/r (r >= R) is now derived in gr.schwarzschild-line-element from DA09 eq. 1 at a = 0. The PPP15 misprint remains in `meaning_source`. |
| `\cos\psi` | add Bel02 eq. 11 use: mu = cos psi with psi the total swept angle of an image ray, which may exceed pi (gr.observed-flux-element-abs-jacobian). |

## Notation bridges (gr_geodesic-r1 cross-entry bullet)

Every source spelling below appears in `truth/entries/gr_geodesic.toml` `latex` fields, which keep
the paper's own symbols. Read them through this table. Each bridge is also recorded in the
affected entries' `notes` (prefixed "Notation bridge (gr_geodesic-r1 cross-entry)").

| source spelling | where | means here | project spelling | collides with (keeps the bare symbol) |
|---|---|---|---|---|
| `E` | DA09 eqs. 6-10; all gr.timelike-* | conserved specific energy of a timelike geodesic, -u_t; in CGS the conserved energy divided by the rest energy m c^2, dimensionless | `\tilde E` | photon energy `E` = hbar omega (rt.notation-val) |
| `R`, `R(r)` | DA09 eq. 6; gr.timelike-radial-equation-schwarzschild, gr.timelike-proper-time-quadrature, gr.timelike-periastron-cubic | radial function of the separated geodesic equations, r^4 (dr/dtau)^2 in geometric units | `\mathcal{R}(r)` | stellar (circumferential) radius `R`, which keeps its meaning in gr.bending-integral-finite-radius and gr.observed-flux-element-abs-jacobian |
| `\rho`, `\rho^2` | DA09 eq. 2; gr.kerr-line-element-bl | Kerr function r^2 + a^2 cos^2 theta (r^2 at a = 0) | `\rho_K^2` | mass density `\rho` |
| `\lambda` | DA09 eqs. 4, 8-10; gr.timelike-* | affine parameter; for delta_1 = 1, lambda = tau + C; geometric lambda_geom = c tau (+ const), a length | `\lambda_{aff}`, and `\tau_{prop}` for proper time | wavelength `\lambda` (colour.notation-cie-cmf) |
| `\theta` | DA09 eqs. 1, 7, 9; gr.kerr-line-element-bl, gr.schwarzschild-line-element, gr.timelike-constants-schwarzschild | Boyer-Lindquist / Schwarzschild polar coordinate (orbit plane theta = pi/2) | `\theta_{Sch}` | RT ray angle `\theta` to the outward normal, mu = cos theta |
| `\phi` | DA09 eqs. 1, 9; gr.timelike-dphi-dtau-schwarzschild, gr.timelike-proper-time-quadrature, gr.timelike-speed-at-infinity | Schwarzschild azimuth in the orbit plane (= `\phi_{orb}` up to origin for the equatorial orbit) | `\phi_{Sch}` | Bel02 surface-element azimuth about D, `\phi` (used as such in gr.observed-flux-element-abs-jacobian) |
| `\mu` | DA09 Sect. 3 (mu = cos theta) | cosine of the Boyer-Lindquist polar coordinate | `\cos\theta_{Sch}` (latex only, no row) | RT direction cosine `\mu` |
| `\mu` | Bel02 eq. 11; gr.observed-flux-element-abs-jacobian | cos psi of the total swept angle of an image ray (psi may exceed pi) | `\cos\psi` | RT direction cosine `\mu` |
| `\psi` | Bel02 eq. 9 with finite upper limit; gr.bending-integral-finite-radius | finite camera-referenced swept angle | `\psi_{cam}` | dictionary `\psi` (infinity-referenced) |
| `M` | DA09 (M = 1), Pech83 | geometric mass M = G M_phys/c^2, a length (= r_g/2) | `GM/c^2` | physical mass `M` in grams |
| `L`, `L_z` | DA09 eqs. 5-9; gr.timelike-* | specific angular momentum about the orbit-plane normal; L_geom = L_CGS/c; \|L_z\| is the total magnitude; L and L_z are signed synonyms only after choosing the equatorial plane | `\tilde L` | bare `L` (retired) |

Unit conversions kept explicit in the entries: M_geom = G M_phys/c^2, lambda_geom = c tau,
L_geom = L_CGS/c, 2M/r -> r_g/r.

Row count: 17 new rows, 5 changed rows.
