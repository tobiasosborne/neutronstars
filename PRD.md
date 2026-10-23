# PRD — Neutron Star Flyby

**Status:** draft v0.1, 2026-10-23. Owner: Tobias Osborne.

## Goal

A video and stills of a free-fall flyby of a real neutron star, RX J1856.5-3754,
against the real galactic sky, that are indistinguishable from what a camera on
that trajectory would record, to the limit of current knowledge. Every pixel is
radiance computed from exact general relativity, measured sky data, and
validated surface physics. Where knowledge runs out, the uncertainty is shown,
not hidden.

## The shot

- **Camera:** on a hyperbolic timelike geodesic of the Schwarzschild spacetime of
  RX J1856. Periastron distance is a parameter, default about 100 km. The camera
  starts and ends far enough out that the star is an unresolved point.
- **Frame:** the camera's local tetrad. Aberration, Doppler, and gravitational
  shifts of the star and the sky follow from that one transformation.
- **Pointing:** parameter, either tracking the star or fixed inertial.
- **Time:** rendered in camera proper time. The close passage lasts milliseconds,
  so a documented time remapping is applied for display and stated on screen.
- **Primary output:** HDR spectral radiance per pixel in absolute units, with
  Stokes parameters. Display images and video are derived with a stated camera,
  exposure, glare, and tone model. Realism lives in the radiance, not the LUT.

## Scope

**In:** exact Schwarzschild geodesics from a finite-radius moving camera;
surface emission as spectral, angular, and polarised radiance from (a)
blackbody, (b) hydrogen atmosphere, (c) thin hydrogen layer over a condensed
surface; surface temperature and field maps constrained by the observed pulse
profile; scale-adaptive star rendering (point with PSF to resolved disc);
background sky from Gaia DR3 re-projected to the star's location plus a
starless diffuse Milky Way; explicit exposure and veiling-glare model.

**Out (for now):** planets, accretion, magnetospheric emission, Kerr or oblate
geometry (the star spins once per 7 s), interstellar extinction at the camera,
anything decorative that is not in the data.

## Requirements

| # | Requirement | Acceptance |
|---|-------------|------------|
| R1 | Exact null geodesics and tetrads | Flat-space limit exact; deflection and visible fraction match Pechenick+1983 and Beloborodov 2002 to stated tolerance |
| R2 | Relativistic camera | Aberration and Doppler reduce to special relativity far from the star; invariance of I_ν/ν³ verified numerically |
| R3 | Real sky | Re-projected Gaia positions and colours match the catalogue within stated tolerance; diffuse map provenance recorded |
| R4 | Star radiance, three surface models | Each model passes an analytic limit and a published-figure comparison; models are interchangeable behind one interface |
| R5 | Earth-distance photometry | A render from 123 pc reproduces the HST optical/UV magnitudes of RX J1856 within published errors, for at least one surface model |
| R6 | HDR output and display model | Absolute units round-trip; exposure, glare, and tone transforms documented and reversible |
| R7 | Provenance | Every equation and constant carries a byte-identical quote from `refs/` |
| R8 | Deliverables | One 4k still at periastron, one video of at least 60 s, and the HDR radiance behind both |

## Success criteria

The Earth-distance render passes R5. A reviewer with the refs can trace every
number in the periastron frame to a source. The periastron still shows the
lensed star ring, limb darkening, and the pole-to-equator gradient with no
element that is not in the data. The surface-state uncertainty is presented
as an ensemble of labelled renders, not a single picture.

## Milestones

1. Spacetime, camera, blackbody disc, real sky. First honest image.
2. Earth-distance photometry harness. Blackbody will fail R5; that failure
   defines the surface-model work.
3. Surface models: audit and salvage the v1 atmosphere solver; add the
   condensed-surface model; pass R5 with at least one model.
4. Flyby video with time remapping, exposure policy, and glare.
5. Uncertainty ensemble and release.

## Open decisions

- Language. Recommendation: Julia, to salvage the v1 microphysics after audit.
- Periastron distance and pointing rule for the hero render.
- Which starless diffuse sky map to use, and its licence.

## Risks

Condensed-surface emissivity needs sources not yet in `refs/` (van Adelsberg
2005, Potekhin 2012, Ho 2007). The 15-decade dynamic range makes the display
model a scientific choice. Gaia re-projection and spectral colours are a data
engineering task that must not swallow the physics time.
