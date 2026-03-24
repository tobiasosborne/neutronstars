# CLAUDE.md — Project Neutron Star

## Mission

Build an end-to-end, open-source, fully reproducible pipeline that computes a **spectral image cube** of a neutron star: the full electromagnetic spectrum (radio through hard X-ray) at every pixel of a resolved image, built from first-principles physics with zero black-box table dependencies. Every equation traces to a locally-stored published paper. Every approximation is named and recorded. Every intermediate quantity can be regenerated from fundamental physics.

The final data product is `I(x, y, ν)` — specific intensity as a function of image-plane position and frequency — from which any visualisation (human vision, false-colour X-ray, polarisation map, single-band, hardness ratio) is a post-processing operation.

Target: Julia-native. Compromises accepted for the tracer bullet (Fortran via `ccall`, downloaded C codes for verification), but the production pipeline must be Julia from the Schrödinger equation to the final pixel.

---

## Sacred Principles

### 1. GROUND TRUTH = LOCAL STRING MATCH

Every equation used in the simulation must match a specific equation in a locally-downloaded published paper or book stored in `refs/`. The match (paper, equation number or page, surrounding context) is recorded in the semantic graph via `af`. No exceptions. "Standard" or "well-known" is not a citation.

Acceptable sources: peer-reviewed papers (arXiv+journal), textbooks with ISBN, institutional theses. Not acceptable: Wikipedia, lecture notes, blog posts, undocumented code, tables from personal websites used as production inputs.

### 2. FAIL FAST, FAIL LOUD

Every function validates inputs and screams on failure. Silent NaN propagation is forbidden. Every numerical integration reports estimated error. Every iterative solver reports convergence status. Every interpolation flags extrapolation. Build assertions into the computation, not around it.

```julia
function opacity_ff_xmode(ν, T, ρ, B)
    @assert T > 0 "T must be positive, got $T"
    @assert isfinite(ν) && ν > 0 "ν must be positive finite, got $ν"
    result = _compute(ν, T, ρ, B)
    @assert isfinite(result) "opacity_ff_xmode returned $result for ν=$ν T=$T ρ=$ρ B=$B"
    return result
end
```

### 3. SKEPTICISM = SUCCESS

Never trust your own previous test results without re-running. Never trust a subagent's claim without independent verification. Never trust a library's output without a sanity check against a known analytic limit. Never trust a downloaded file without checksumming it.

Every physics module requires at minimum TWO independent checks: one analytic limit, one numerical comparison against a published value with the source cited. Record all verification in `verification/`.

### 4. BE BOLD, LEARN FROM FAILURE

Failing verbosely and quickly is success. Every failure must produce a clear error message, a hypothesis about why, and a note in `notes/failures.md`. Getting a crude end-to-end tracer bullet working beats a perfect submodule with nothing to plug into.

### 5. ZERO EXTERNAL TABLE DEPENDENCIES

Every numerical quantity in the production pipeline is computed from first principles. Published tables are used ONLY for verification (checking your code reproduces known results), NEVER as production inputs. The pipeline must be fully reproducible from: the Schrödinger equation, the Dirac equation, Maxwell's equations, statistical mechanics, and general relativity.

This means recomputing from scratch:
- OCP Coulomb free energy (classical Monte Carlo, ~150 lines)
- Magnetic hydrogen energy levels (variational diagonalisation, ~250 lines)
- Oscillator strengths (dipole matrix elements, ~50 lines on top of the above)
- Photoionisation cross-sections (continuum wavefunctions + matrix elements, ~100 lines)
- Free-free Gaunt factors (Coulomb wavefunctions, or cross-check against van Hoof's GPL code)
- Magnetic free-free opacities (from published formulae, not from Potekhin's Fortran)
- Non-magnetic and magnetic atmosphere spectra (Feautrier RT solver)

Existing codes (Potekhin Fortran routines, McPHAC, HFFERII, van Hoof's Gaunt code) are verification targets, not dependencies.

---

## Output Data Product: Spectral Image Cube

The primary output is NOT an RGB image. It is a spectral image cube.

### Specification

```julia
struct SpectralImageCube
    nx::Int; ny::Int; nν::Int
    ν_grid::Vector{Float64}           # observer-frame frequencies [Hz], log-spaced
    pixel_scale::Float64              # angular size per pixel [rad]
    I::Array{Float64, 3}             # (nx, ny, nν) specific intensity [erg/s/cm²/Hz/sr]
    Q::Union{Nothing, Array{Float64,3}}  # Stokes Q (optional)
    U::Union{Nothing, Array{Float64,3}}  # Stokes U (optional)
    # Per-pixel diagnostic metadata
    hit_fraction::Matrix{Float64}     # fraction of ray bundle hitting surface
    mean_redshift::Matrix{Float64}
    mean_T::Matrix{Float64}
    mean_B::Matrix{Float64}
    # Provenance
    params::NSParameters
    n_rays_per_pixel::Int
    atmosphere_id::String
    code_version::String
    approximations::Vector{String}    # active approximation node IDs from af graph
end
```

### Frequency grid

~500–1000 points logarithmically spaced from 10¹³ Hz (far IR) to 10²⁰ Hz (hard X-ray). Log spacing naturally concentrates resolution in the X-ray band (10¹⁶–10¹⁸ Hz) where thermal emission peaks and atmospheric features appear. The optical band (4–7.5 × 10¹⁴ Hz) gets sparse but sufficient sampling since the spectrum is smooth Rayleigh-Jeans there.

### Per-pixel ray bundle

Each pixel gets M rays (default M=16) sampling its solid angle. For each ray:

1. **GR geodesic** → hit/miss, surface coordinates (θ,φ), emission angle cos θ_e, redshift (1+z), Doppler factor δ
2. **Surface model** → local T, local B vector, angle between B and surface normal
3. **Atmosphere spectrum** → I_ν^(X)(ν'), I_ν^(O)(ν') in both polarisation modes at emission-frame frequencies
4. **(Optional) Magnetosphere RT** → modify spectrum along the ray path through extended plasma
5. **Relativistic transform** → observer-frame I_ν^obs(ν_obs) = I_ν'(ν') / ((1+z)δ)³ with ν' = (1+z)δ ν_obs
6. **Accumulate** into pixel with solid-angle weight, sum polarisation modes to get Stokes I (keep separate for Q, U)

Rays that miss the star contribute zero.

### Rendering (post-processing on the cube)

```julia
render_human_vision(cube)              # CIE 1931 XYZ → sRGB → tone map
render_false_colour(cube, bands)       # arbitrary band → RGB mapping
render_bolometric(cube)                # ∫ I_ν dν per pixel
render_mono(cube, ν)                   # monochromatic slice
render_spectrum(cube, x, y)            # full SED at one pixel
render_band_ratio(cube, band1, band2)  # hardness ratio map
render_polarisation(cube)              # polarisation degree and angle
```

---

## Physics Pipeline: Six Modules

### Module 1: Neutron Star Structure (TOV + EOS)

**Input:** central density ρ_c, equation of state.
**Output:** M, R, g_surface, metric functions.

Solve the Tolman-Oppenheimer-Volkoff equations. For the EOS, implement the Potekhin et al. (2013) analytical fits to the BSk24 unified EOS — these are published formulae, not tables. The fit coefficients are printed in the paper.

**Key references:**
- Shapiro & Teukolsky (1983) Ch 5-6: TOV equations
- Potekhin et al., A&A 560, A48 (2013): analytical EOS fits with explicit coefficients

**Verification:** reproduce Potekhin+(2013) Table 1 (M, R for given EOS) to 4 significant figures. Newtonian limit for M→0.

### Module 2: Opacities κ_ν^(j)(ν, T, ρ, B, composition)

The hardest module. Three sub-components:

**2a. Free-free (bremsstrahlung) absorption — magnetic case.**
Both X-mode and O-mode cross-sections from Potekhin & Chabrier (2003) Appendix A. These are explicit formulae involving thermal averages of Coulomb matrix elements between Landau-level states. In the strongly quantising regime (ℏω_c ≫ kT), most electrons occupy the ground Landau level and the expressions simplify.

For the B=0 limit: Kramers formula with Gaunt factor corrections. Recompute Gaunt factors from Coulomb wavefunctions (or verify against van Hoof et al. 2014 GPL code at `data.nublado.org/gauntff/`).

**Key references:**
- Potekhin & Chabrier, A&A 399, 1007 (2003): magnetic ff opacity, Appendix A equations
- Rybicki & Lightman (1979) Ch 5: non-magnetic bremsstrahlung
- van Hoof et al., MNRAS 444, 420 (2014): Gaunt factors from first principles, source code GPL

**2b. Scattering — magnetic Thomson/Compton.**
Analytic QED cross-sections from Mészáros (1992). Polarisation-dependent, anisotropic. X-mode is suppressed by (ν/ν_c)⁴ below cyclotron frequency. Compton correction via Kompaneets equation.

**Key references:**
- Mészáros, "High-Energy Radiation from Magnetized Neutron Stars" (1992) Ch 3-5
- Rybicki & Lightman (1979) Ch 7: Kompaneets equation

**2c. Bound-free (photoionisation) — only if partially ionised.**
Requires magnetic hydrogen atomic data. This is the dependency chain:

*Energy levels:* Solve the Schrödinger equation for H in uniform B-field. Adiabatic approximation: transverse wavefunctions are Landau functions (analytic), longitudinal problem is a 1D Schrödinger equation with effective Coulomb potential. For higher accuracy: expand in Landau basis × longitudinal functions, diagonalise the Hamiltonian matrix. Multiple independent codes exist for verification (HFFERII in CPC archive, Cao et al. 2019 in CPC archive, Xi et al. 1992, Liu et al. 2006). Verify against Ruder et al. (1994) Table 5.1 and Schimeczek & Wunner (2014) published data.

*Oscillator strengths:* Dipole matrix elements ⟨ψ_f|r⃗|ψ_i⟩ between the wavefunctions computed above. Factorises into transverse (analytic for Landau functions) × longitudinal (1D numerical quadrature) parts.

*Photoionisation cross-sections:* Same matrix element but with the final state being a continuum wavefunction (ODE solve for the 1D Coulomb scattering state at E>0). Source: Potekhin & Pavlov (1997), Potekhin et al. (2004).

**Key references:**
- Ruder et al., "Atoms in Strong Magnetic Fields" (1994) ISBN 3-540-57499-0
- Schimeczek & Wunner, Comp. Phys. Comm. 185, 2992 (2014): HFFERII code description
- Cao et al., Comp. Phys. Comm. 241, 129 (2019): independent H-in-B solver
- Potekhin & Pavlov, ApJ 483, 414 (1997): photoionisation cross-sections

**Decision point:** Start with fully ionised hydrogen (T_eff > 5×10⁵ K). This eliminates the bound-free dependency entirely. Partial ionisation is Phase 4.

### Module 3: Atmosphere Radiative Transfer

**Input:** T_eff, log g, B, composition, frequency grid.
**Output:** I_ν^(j)(ν, cos θ_e) — emergent specific intensity in each polarisation mode as a function of frequency and emission angle.

Solve two coupled radiative transfer equations (one per polarisation mode) in a plane-parallel atmosphere in hydrostatic and radiative equilibrium. Numerical method: Feautrier scheme (Mihalas 1978 Ch 6). Iterate temperature profile to enforce constant integrated flux (temperature correction / ALI).

**Implementation order:**
1. Non-magnetic (B=0) pure H, fully ionised. Verify against McPHAC (<1% across 0.01–10 keV).
2. Add magnetic opacities (two modes). Verify B→0 limit recovers step 1.
3. Verify magnetic spectra against Suleimanov et al. (2009) published figures.
4. (Later) Add partial ionisation using Module 2c.

**Key references:**
- Mihalas, "Stellar Atmospheres" (1978): Feautrier method, boundary conditions, temperature correction
- Haakonsen et al., ApJ 749, 52 (2012): McPHAC paper (verification target)
- McPHAC source code: github.com/McPHAC/McPHAC (structural reference)
- Suleimanov, Potekhin & Werner, A&A 500, 891 (2009): magnetic atmosphere scheme + verification spectra

**Pre-computation:** Once the solver works, generate a grid of atmosphere models over (T_eff, log g, B, cos θ_e) and store as your own reproducible table (HDF5 with full provenance metadata). The ray tracer interpolates from this grid at runtime.

### Module 4: Surface Model

**Input:** NS parameters (M, R, B_pole, magnetic obliquity α, T_core or age).
**Output:** T(θ,φ), B⃗(θ,φ) across the stellar surface.

**Temperature map:** Greenstein & Hartke (1983) analytic model for a dipole magnetic field:
T(θ_B) = T_eq + (T_pole − T_eq) cos²θ_B

For more realism: use NSCool or MATINS-derived envelope models (Potekhin, Pons & Page 2015) relating internal temperature to surface temperature as a function of local B. But start with Greenstein-Hartke.

**Magnetic field geometry:** Dipole:
B_r = B_pole cos θ_B, B_θ = (B_pole/2) sin θ_B, |B| = (B_pole/2)√(1+3cos²θ_B)

where θ_B is the magnetic colatitude (accounting for obliquity α between magnetic and rotation axes).

**Key references:**
- Greenstein & Hartke, ApJ 271, 283 (1983)
- Potekhin, Pons & Page, SSRv 191, 239 (2015) / arXiv:1507.06186

### Module 5: GR Ray Tracing (Geodesics)

**Input:** observer position, image-plane resolution, NS spacetime (M, R, spin if Kerr).
**Output:** per ray — hit/miss, surface coordinates, emission angle, redshift, Doppler factor, path through exterior spacetime.

**Schwarzschild (non-rotating):** The geodesic equation is integrable via elliptic integrals. Map impact parameter b to deflection angle, invert to get surface coordinates and emission angle. This is fast — no ODE integration needed.

**Kerr (rotating):** Dexter & Agol (2009). Save for Phase 4.

The key outputs for each ray are:
- Whether it hits the surface (impact parameter < critical value)
- Surface hit coordinates (θ, φ) in the NS frame
- Local emission angle cos θ_e (angle between ray and surface normal, in the local frame)
- Total gravitational redshift 1+z = (1−2GM/Rc²)^{−1/2}
- Doppler boost δ = 1/[γ(1−β⃗·n̂)] if rotating
- Path segments through exterior spacetime (for magnetosphere RT)

**Key references:**
- Pechenick, Ftaclas & Cohen, ApJ 274, 846 (1983): exact image of NS in Schwarzschild
- Beloborodov, ApJ 566, L85 (2002): simplified + exact formulae
- Dexter & Agol, ApJ 696, 1616 (2009): Kerr geodesics

**Verification:**
- Flat space limit (M→0): visible fraction → exactly 1/2
- Reproduce Pechenick+(1983) Fig 1
- Cross-check against X-PSI's geodesic module (github.com/xpsi-group/xpsi, Python/C, open source)
- Total observed luminosity consistency: L_obs = L_emit / (1+z)²

### Module 6: Colorimetry and Rendering

**Input:** spectral image cube I(x,y,ν).
**Output:** RGB images under various rendering modes.

**Human vision rendering:**
1. Integrate I_ν against CIE 1931 2° observer colour matching functions x̄(λ), ȳ(λ), z̄(λ). Source: CIE 015:2004, data from `cvrl.ioo.ucl.ac.uk`.
2. XYZ → linear sRGB via the IEC 61966-2-1 matrix (exact published values).
3. sRGB gamma: C_srgb = 12.92C if C ≤ 0.0031308, else 1.055C^{1/2.4} − 0.055.
4. Tone mapping: Reinhard global operator or ACES filmic. The dynamic range is enormous (poles vs equator differ by orders of magnitude in optical).

**False colour:** Map arbitrary frequency bands to RGB channels. User-configurable.

**Verification:**
- 5778K blackbody → solar white in sRGB
- CIE standard illuminant D65 chromaticity coordinates
- XYZ→sRGB matrix matches IEC 61966 to machine precision

**Key references:**
- CIE 015:2004 (colour matching functions)
- IEC 61966-2-1:1999 (sRGB specification)
- Reinhard et al., SIGGRAPH 2002 (tone mapping)

---

## Equation of State: The Coulomb Correction

The EOS for the atmosphere requires the free energy of a magnetised, partially degenerate electron-ion plasma. The ideal gas + Landau quantisation parts are analytic (Fermi-Dirac integrals, published formulae in Potekhin & Chabrier 2001). The Coulomb interaction correction depends on Monte Carlo data for the classical one-component plasma (OCP).

**Recompute the OCP from scratch:**
- Classical MC with Ewald summation + Metropolis sampling. N ≥ 1000 particles, scan Γ = 1 to 300.
- Thermodynamic integration to get free energy from internal energy.
- Verify against Slattery, Doolen & DeWitt (1982) Table I. Also verify against modern recalculations: Caillol (1999), Kozhberov (2025, arXiv:2511.04516).

**Key references:**
- Slattery, Doolen & DeWitt, Phys. Rev. A 26, 2255 (1982): original OCP MC
- Potekhin & Chabrier, Phys. Rev. E 62, 8554 (2000): fit to OCP data used in NS EOS
- Frenkel & Smit, "Understanding Molecular Simulation" Ch 12: Ewald summation recipe
- Kozhberov, arXiv:2511.04516 (2025): modern independent MD recalculation

This is ~150 lines of Julia. The original code ran on a CDC 7600 with 490 KB of memory. A modern laptop is embarrassingly overpowered.

---

## Existing Code Landscape (verification targets, not dependencies)

No monolith project covers this pipeline. The field is fragmented:

| Code | What it does | What it DOESN'T do | Language | Availability |
|------|-------------|-------------------|----------|-------------|
| X-PSI | GR ray tracing, pulse profiles, Bayesian inference | Compute atmosphere spectra (loads FITS tables) | Python/C | GitHub, open |
| McPHAC | Non-magnetic H atmosphere RT | Magnetic case, ray tracing, rendering | C | GitHub, open |
| NSCool | 1D cooling evolution T(age) | Atmosphere spectra, ray tracing, rendering | Fortran 77 | Web, free |
| MATINS | 3D magnetothermal evolution, surface T maps | Atmosphere spectra, ray tracing, rendering | Fortran/C++ | NOT public |
| HFFERII | Magnetic H atom energy levels + oscillator strengths | Opacities, RT, ray tracing | Fortran | CPC archive |
| Cao et al. (2019) | Independent magnetic H atom solver | Same | Fortran | CPC archive |
| van Hoof | B=0 Gaunt factors from Coulomb wavefunctions | Magnetic case | C/Fortran | GPL, web |
| Potekhin routines | Microphysics: conductivity, EOS, opacity fits | Everything else | Fortran | ioffe.ru (fragile) |
| NSA/NSATMOS/NSMAX | Interpolate pre-computed atmosphere FITS tables | Generate the tables | C (XSPEC) | HEASARC |
| Suleimanov code | Magnetic atmosphere RT (modified ATLAS) | Everything else | Fortran | NOT public |
| CompOSE | Database of EOS tables | Compute them (just stores them) | Web DB | compose.obspm.fr |

**The two things nobody has built:**
1. An open-source magnetic atmosphere RT solver
2. Visual/spectral rendering from a spectral image cube

These are the core contributions of this project.

---

## Semantic Graph (af tool)

Use `af` to build a directed graph of claims, equations, and approximations. Every node links to a local PDF in `refs/`.

**Node types:**

- **EQUATION**: A mathematical relation. Fields: LaTeX, source paper, equation number, page, verified (bool), implemented_in (file:line).
- **APPROXIMATION**: A simplification. Fields: description, validity range, source, affects (list of equation nodes), impact assessment.
- **CLAIM**: A physics assertion. Fields: statement, source, status (accepted/contested/placeholder).

**Required edges:** DEPENDS_ON, APPROXIMATES, VERIFIED_BY, IMPLEMENTS, DERIVED_FROM.

**Minimum graph size before writing simulation code:** ≥50 equation nodes, ≥20 approximation nodes, all with source links to local PDFs.

---

## Phased Development

### Phase 1: Literature + Semantic Graph

Download all references to `refs/`. Build `refs/INDEX.md` with SHA256 checksums. Construct the `af` semantic graph. Produce `notes/approximations.md` (master list of all approximations with validity ranges) and `notes/decisions.md`.

**Key decision to record:** "Restrict to fully ionised hydrogen atmosphere for Phases 2-3. Valid for T_eff > 5×10⁵ K. Eliminates bound-free opacity and magnetic H atom dependency. Partial ionisation deferred to Phase 4."

### Phase 2: Tracer Bullet (end-to-end with simplified atmosphere)

Get photons from the NS surface to a spectral image cube using a modified blackbody placeholder for the atmosphere: I_ν = f_col⁴ B_ν(f_col T, ν) × limb_darkening(cos θ_e). Record this as APPROXIMATION node "colour_corrected_blackbody" with flag PLACEHOLDER.

Deliverables:
- Working TOV solver (verified against Potekhin+ 2013)
- Working Schwarzschild ray tracer (verified against Pechenick+ 1983)
- Working surface model (Greenstein-Hartke)
- Working spectral image cube with full frequency grid
- Working rendering: human vision (dim blue-white disk) + false-colour X-ray (revealing temperature structure + lensing)
- All verification checks passing with printed comparison values

### Phase 3: Real Atmosphere

Replace the blackbody placeholder with a self-consistent Feautrier RT solver.

**3a.** Non-magnetic atmosphere. Verify against McPHAC to <1%.
**3b.** Magnetic opacities (ff + scattering, both modes). Verify B→0 limit.
**3c.** Full magnetic atmosphere. Verify against Suleimanov+(2009) figures.
**3d.** Pre-compute atmosphere grid. Store as HDF5 with full provenance.
**3e.** Re-render spectral image cube with real atmosphere.

### Phase 4: Extensions

In priority order:
1. OCP Monte Carlo (replace any residual fits with your own simulation data)
2. Partial ionisation (magnetic H atom solver + bound-free opacities)
3. Vacuum polarisation / QED mode conversion (Ho & Lai 2001; matters for B > 10¹³ G)
4. Kerr ray tracing (Dexter & Agol 2009; adds rotation, Doppler, frame dragging)
5. Magnetosphere volume rendering (Lyr.jl MC engine; pair plasma, cyclotron features)
6. Condensed surface emission (van Adelsberg et al. 2005; alternative to gaseous atmosphere)
7. ISM absorption (Wilms, Allen & McCray 2000; for specific observed NS at known distance)

---

## Verification Matrix

Every module must pass ALL checks before being marked complete.

| Module | Analytic Limit | Literature Value | Cross-Code |
|--------|---------------|-----------------|------------|
| TOV | Newtonian limit (M→0) | Potekhin+(2013) Table 1 | NSCool |
| EOS (B=0) | Ideal gas (T→∞) | Potekhin & Chabrier (2001) | Potekhin Fortran |
| EOS (B>0) | B→0 recovers non-mag | Potekhin & Chabrier (2001) | Potekhin Fortran |
| OCP MC | Debye-Hückel (Γ→0) | Slattery+(1982) Table I | Kozhberov (2025) |
| Gaunt factors | Kramers (low ν) | van Hoof+(2014) tables | van Hoof GPL code |
| Opacity ff (B=0) | Kramers at low ν | Rybicki & Lightman examples | McPHAC |
| Opacity ff (B>0) | B→0 recovers non-mag | Potekhin & Chabrier (2003) Fig 1 | Potekhin Fortran |
| Opacity scat (B=0) | σ = σ_T | σ_T = 6.6524e-25 cm² | McPHAC |
| Opacity scat (B>0) | ν≪ν_c: X-mode suppressed | Mészáros (1992) Ch 5 | — |
| Atmosphere (B=0) | Grey: Hopf function | McPHAC spectra (<1%) | McPHAC |
| Atmosphere (B>0) | B→0 recovers non-mag | Suleimanov+(2009) Figs 4-6 | XSPEC nsmax tables |
| Geodesics | M→0: flat space | Pechenick+(1983) Fig 1 | X-PSI |
| Surface model | Uniform T,B: isotropic | Greenstein & Hartke (1983) | — |
| CIE colorimetry | 5778K BB → solar white | CIE standard illuminants | — |
| sRGB transform | Identity for in-gamut | IEC 61966 spec | — |
| Flux conservation | L_obs = L_emit/(1+z)² | Analytic | X-PSI |

---

## Coding Standards

- **CGS throughout** the physics code. Define `PhysicalConstants` module with CODATA 2018 values, each citing the source.
- **Every physics function** gets a docstring citing its source: equation number, paper, page.
- **No magic numbers.** Every constant gets a name and a citation.
- **~200 lines per file.** Decompose aggressively.
- **Test files mirror source:** `src/opacity/magnetic_ff.jl` → `test/opacity/test_magnetic_ff.jl`
- **Verification prints comparison values:** expected (with source), computed, relative error, tolerance. "Test passed" alone is not acceptable.

---

## Anti-Patterns

1. **Implementing an equation from memory** without finding it in a local PDF.
2. **Downloading a table and using it as a production input** rather than recomputing from published formulae.
3. **"Test passed"** without printing expected value, computed value, error, and source.
4. **Silent fallbacks.** Clamping, extrapolating, or substituting defaults without logging.
5. **"It looks right."** Visual inspection is not verification. Quantitative comparison against a published number is.
6. **Skipping the semantic graph.** If you can't trace a number back to a local PDF through the graph, the number is not trustworthy.
7. **Premature optimisation.** Get it correct first, then fast. Profile before optimising.

---

## Directory Structure

```
neutron-star/
├── CLAUDE.md
├── refs/                    # Downloaded papers, books, code
│   ├── INDEX.md             # bibkey → filename → what it provides → SHA256
│   ├── code/                # Reference codebases (McPHAC, van Hoof, etc.)
│   └── *.pdf
├── graph/                   # af semantic graph
├── src/
│   ├── NeutronStar.jl       # Top-level module
│   ├── constants.jl         # Physical constants with CODATA citations
│   ├── structure/           # TOV solver, EOS
│   ├── atomic/              # Magnetic H atom, energy levels, wavefunctions
│   ├── plasma/              # OCP Monte Carlo, Coulomb corrections
│   ├── opacity/             # ff, scattering, bf — magnetic and non-magnetic
│   ├── atmosphere/          # Feautrier RT solver, temperature iteration
│   ├── surface/             # T(θ,φ), B(θ,φ) models
│   ├── geodesics/           # Schwarzschild (and later Kerr) ray tracing
│   ├── magnetosphere/       # (Phase 4+) Volume RT
│   ├── cube/                # SpectralImageCube type and operations
│   ├── render/              # Colorimetry, false colour, tone mapping
│   └── pipeline/            # End-to-end orchestration
├── verification/            # Quantitative checks against literature
│   └── VERIFICATION_LOG.md
├── notes/
│   ├── failures.md
│   ├── decisions.md
│   └── approximations.md
├── output/                  # Spectral cubes and rendered images
├── Project.toml
└── test/                    # Mirrors src/ structure
```

---

## First Session Checklist

1. Create directory structure.
2. Begin downloading arXiv papers: `wget https://arxiv.org/pdf/XXXX.XXXXX -O refs/author_year.pdf`
3. Clone McPHAC: `git clone https://github.com/McPHAC/McPHAC.git refs/code/McPHAC`
4. Clone X-PSI: `git clone https://github.com/xpsi-group/xpsi.git refs/code/xpsi`
5. Download van Hoof Gaunt factors: `wget https://data.nublado.org/gauntff/gauntff.dat -O refs/data/gauntff.dat`
6. Initialise Julia project.
7. Implement `src/constants.jl` with CODATA citations.
8. Begin `af` graph construction with top-level claims.
9. Implement TOV solver (simplest self-contained module with clear verification).
10. **Do not start atmosphere code until Phase 1 graph has ≥50 equation nodes.**
