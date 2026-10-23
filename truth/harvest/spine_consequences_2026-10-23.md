# Spine consequences from the 2026-10-23 Milestone-1 harvest (compiled from lane reports)

Entry status as of writing: see `python3 -I truth/check.py`. Ids below exist (blessed or proposed) unless marked.

## pipeline.constant-set (blessed; extend depends_on)
add: const.au, const.parsec, const.gm-sun-nominal, const.gm-sun-nominal-cgs, const.m-sun-nominal-cgs, const.erg-joule, const.gauss-tesla.
(const.esu-from-coulomb is REJECTED r2; const.r-sun-nominal, const.l-sun-nominal dropped.) Statement should say the Gaussian bridge is a project definition pending Jackson.

## pipeline.ns-exterior-spacetime (blessed; extend)
depends_on = ["gr.schwarzschild-line-element", "gr.schwarzschild-radius", "gr.approx-schwarzschild-exterior", "gr.apparent-quantities", "surface.ppp15-redshifted-quantities", "surface.phc16-redshift-relations", "pipeline.constant-set"]
Also in gr.toml: gr.approx-schwarzschild-exterior and gr.frequency-redshift currently depend on the REJECTED gr.metric-static-spherical; repoint to gr.schwarzschild-line-element (blessed).

## pipeline.camera-worldline (REJECTED r1; re-propose)
depends_on = ["gr.schwarzschild-line-element", "gr.timelike-constants-schwarzschild", "gr.timelike-radial-equation-schwarzschild", "gr.timelike-dt-dtau-schwarzschild", "gr.timelike-dphi-dtau-schwarzschild", "gr.timelike-proper-time-quadrature", "gr.timelike-periastron-cubic", "gr.timelike-angular-momentum-from-periastron", "gr.timelike-unbound-condition", "gr.timelike-speed-at-infinity", "pipeline.ns-exterior-spacetime"]
Drop gr.static-observer-four-velocity. Statement must fix: orientation convention (sign of L_z, orbit plane = equatorial theta_Sch = pi/2), start/end radius rule (star unresolved: angular diameter below the PSF first zero, see rt.airy-first-zero and pipeline.star-scale-adaptive), M in grams via const.gm-sun-nominal-cgs and const.m-sun-nominal-cgs, proper-time parametrisation tau_prop with the unsigned two-leg quadrature.

## pipeline.camera-tetrad-projection (REJECTED r1; re-propose)
depends_on = ["gr.static-tetrad", "gr.camera-four-velocity-static-frame", "gr.tetrad-lorentz-transformation", "gr.tetrad-vector-components", "gr.camera-aberration", "gr.camera-projection-convention", "gr.null-geodesic-constants", "gr.observer-sky-angle", "gr.disk-angular-radius", "pipeline.camera-worldline", "pipeline.ns-exterior-spacetime"]
Drop gr.emission-angle-definition. Direction convention: N points toward the apparent source; p = (h nu_cam/c)(e_(0) - N^a e_(a)). Projection: pinhole or equidistant per gr.camera-projection-convention (open PRD decision noted there).

## pipeline.frequency-shift (REJECTED? check status; re-propose or extend)
depends_on = ["gr.frequency-ratio-general", "gr.intensity-g-cubed", "gr.invariant-inu-over-nu3", "gr.frequency-redshift", "gr.doppler-factor-sr", "gr.camera-blueshift-celestial", "gr.static-observer-four-velocity", "gr.disk-uniform-brightness", "pipeline.camera-tetrad-projection", "pipeline.ns-exterior-spacetime", "pipeline.null-geodesic-trace"]

## pipeline.null-geodesic-trace (check status; extend)
depends_on = ["gr.schwarzschild-line-element", "gr.null-geodesic-constants", "gr.null-geodesic-radial", "gr.null-geodesic-components", "gr.orbit-equation-u", "gr.bending-integral", "gr.bending-integral-pechenick", "gr.bending-integral-finite-radius", "gr.photon-sphere", "gr.impact-parameter-max", "gr.max-deflection-integral", "gr.emission-angle-definition", "gr.emission-angle-impact-parameter", "gr.approx-opaque-surface", "gr.observed-flux-element-abs-jacobian", "gr.flat-space-half-visible", "gr.full-visibility-radius", "gr.dphi-max-r16", "gr.circle180-radius-r16", "gr.beloborodov-cosine-relation", "gr.beloborodov-cosine-error-bound", "gr.bending-linear-in-u", "gr.visible-fraction-cosine", "pipeline.camera-tetrad-projection", "pipeline.ns-exterior-spacetime", "pipeline.surface-maps-field"]
(ids that do not exist: drop and say so.) gr.observed-flux-element-abs-jacobian is a far-field verification target only.

## pipeline.surface-maps-field (extend)
add: surface.py01-gr-dipole-magnitude, surface.py01-gr-dipole-field-angle, surface.py01-gr-dipole-a-of-x, surface.ws83-gr-dipole-components, surface.field-colatitude-from-geometry (all blessed). Statement: GR dipole primary; flat tan Theta = 1/2 tan theta_B is its a -> 1/2 limit.

## pipeline.surface-radiance-blackbody (extend)
add: rt.planck-function-bnu, rt.planck-integral-stefan-boltzmann (blessed), const.planck, const.boltzmann, const.speed-of-light. (rt.planck-function-blambda proposed, r3 pending.)

## pipeline.surface-parameters (extend)
add: surface.rxj1856-rinf-ho2007-rescaled-123pc (blessed).

## pipeline.display-colour (extend)
add: colour.srgb-stokes-decode-power, colour.k-cd-luminous-efficacy (blessed). Still missing: tristimulus normalisation (CIE 015, paywalled) -> keep node yellow and say so.

## pipeline.display-exposure-tone (extend)
add: colour.display-exposure-convention, colour.reinhard-simple-inverse (blessed).

## pipeline.star-scale-adaptive (extend)
add: rt.psf-aperture-convention (blessed; it depends on rt.airy-intensity-pattern [proposed, r3 pending] and rt.airy-first-zero [blessed]).

## pipeline.sky-radiance (REJECTED r1; re-propose)
add (sky lane, blessed unless noted): obs.rxj1856-position-walter2010, obs.gaia-reference-epoch-j2016, obs.gaia-reference-frame-crf3, obs.gaia-dr3-astrometric-parameters, obs.gaia-standard-model-of-stellar-motion, obs.gaia-perspective-acceleration, obs.gaia-band-flux-integral, obs.gaia-vegamag-definition, obs.gaia-ab-definition, obs.gaia-zero-points-vegamag (proposed r3), obs.gaia-zero-points-ab, obs.gaia-passband-distribution, obs.gaia-xp-ecs-definition, obs.gaia-xp-ecs-spd-sed-units (proposed r3), obs.gaia-xp-sampled-subset, obs.gaia-xp-absolute-calibration (proposed r3), obs.gaia-magnitude-range-and-bright-end.
add (sky2 lane; status per sky2 r1 apply, check): obs.agis-standard-model-direction, obs.agis-astrometric-parameter-set, obs.bl-propagation-direction-nolt, obs.bl-propagation-parallax-nolt, obs.bl-propagation-pm-nolt, obs.bl-propagation-mur-nolt, obs.bl-lighttime-neglect-criterion, obs.rxj1856-epoch-propagation-estimate, obs.hipparcos-optical-realization-icrs; candidates pending the diffuse-map PRD decision: obs.mellinger-background-and-star-content, obs.mellinger-photometric-calibration, obs.mellinger-availability.
Statement must say: Mellinger is NOT starless and has no licence on disk (open PRD decision: mask+inpaint with Gaia positions / another map / Gaia total minus catalogued flux); stars fainter than G = 15 have no sampled XP spectrum; data files still to fetch: Gaia archive rows, passband tables (VizieR), XP sampled spectra, CALSPEC Vega, Hipparcos rows.
