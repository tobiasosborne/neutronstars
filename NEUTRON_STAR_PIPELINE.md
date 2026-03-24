# CLAUDE.md — Project Neutron Star: Physically Exact Visualisation from First Principles

## Mission

Build a complete, open-source, reproducible computational pipeline that renders a physically accurate image of a neutron star as seen by a nearby observer. Every photon in the final image must be traceable through an unbroken chain of equations back to published, locally-stored, peer-reviewed literature. No black-box tables. No "trust me" constants. No flaky websites. Every approximation must be named, justified, and recorded in a semantic graph.

The final deliverable is a Julia-native codebase that takes as input a set of neutron star parameters (M, R, B, T_core, age, viewing geometry) and produces a multi-band image with full spectral information at every pixel, suitable for both false-colour scientific visualisation and human-vision-accurate rendering.

---

## Sacred Principles

These are non-negotiable. Violating any of them is a build-breaking event.

### 1. GROUND TRUTH = LOCAL STRING MATCH

Every equation used in the simulation must be matched to an equation in a locally-downloaded published paper or textbook. "Locally-downloaded" means a PDF or DjVu file on disk in `refs/`. "String match" means you can point to the specific equation number, page, and surrounding text. The match must be recorded in the semantic graph.

**Acceptable ground truth sources:**
- Published peer-reviewed papers (arXiv PDFs count if the paper is published; note the journal reference)
- Textbooks with ISBN
- Doctoral theses from institutional repositories

**NOT acceptable:**
- Wikipedia
- Lecture notes without formal publication
- Blog posts
- "I know this equation" or "this is standard"
- Tables downloaded from someone's personal website without the paper that derived them
- Any equation not matched to a local file

**Workflow:** Before implementing ANY equation, you must:
1. Identify the source paper/book
2. Download it to `refs/`
3. Find the exact equation (by number or by page + position)
4. Record the match in the semantic graph via `af`
5. Only then implement it

### 2. FAIL FAST, FAIL LOUD

Every function must validate its inputs and scream on failure. Silent NaN propagation is a capital offence. Use Julia's type system aggressively. Assertion-heavy code is good code.

```julia
# YES
function rosseland_opacity(T::Float64, ρ::Float64, B::Float64)
    @assert T > 0 "Temperature must be positive, got $T"
    @assert ρ > 0 "Density must be positive, got $ρ"
    @assert B >= 0 "Magnetic field must be non-negative, got $B"
    result = _compute_opacity(T, ρ, B)
    @assert isfinite(result) "Opacity computation returned $result for T=$T, ρ=$ρ, B=$B"
    return result
end

# NO
function rosseland_opacity(T, ρ, B)
    return _compute_opacity(T, ρ, B)  # silent NaN factory
end
```

Every numerical integration must report its estimated error. Every iterative solver must report convergence status. Every interpolation must report whether it extrapolated. Warnings are not errors, but unhandled warnings accumulate into errors. Build a `DiagnosticLog` that every function writes to.

### 3. SKEPTICISM = SUCCESS

Never trust:
- Your own previous test results (re-run them)
- Any subagent's claim that "this works" (verify independently)
- Any library's output without a sanity check against a known limit
- Any downloaded file's integrity (checksum it)
- Any equation you implemented without testing against at least one known value from the source paper

For every physics module, implement at least TWO independent checks:
- **Limit check:** Does the code reproduce the correct answer in a known analytic limit?
- **Literature check:** Does the code reproduce a specific numerical value or figure from a published paper?

Record all verification results in `verification/` with the exact comparison values and their sources.

### 4. BE BOLD, LEARN FROM FAILURE

This is a hard project. Some parts will fail on the first attempt. That is fine. What is not fine is failing silently or failing without learning. Every failure must produce:
- A clear error message explaining what went wrong
- A hypothesis about why
- A recorded note in `notes/failures.md`

Ambition is rewarded. Getting a crude end-to-end tracer bullet working (even with placeholder physics) is worth more than a perfect opacity module with no pipeline to plug into.

---

## Project Structure

```
neutron-star/
├── CLAUDE.md                    # This file
├── refs/                        # ALL downloaded papers and books (PDFs)
│   ├── INDEX.md                 # Paper index: bibkey → filename → what it's used for
│   └── *.pdf
├── graph/                       # Semantic graph (af tool output)
│   ├── claims/                  # Individual claim files
│   ├── approximations/          # Named approximations with justifications
│   └── equations/               # Equation-to-source mappings
├── src/                         # Julia source code
│   ├── NeutronStar.jl           # Main module
│   ├── eos/                     # Equation of state
│   ├── opacity/                 # Opacity calculations
│   ├── atmosphere/              # Radiative transfer solver
│   ├── surface/                 # Surface temperature map
│   ├── geodesics/               # GR ray tracing
│   ├── colorimetry/             # Spectrum → human colour
│   └── pipeline/                # End-to-end orchestration
├── verification/                # Test results against literature
│   ├── eos_checks.jl
│   ├── opacity_checks.jl
│   ├── atmosphere_checks.jl
│   ├── geodesic_checks.jl
│   └── VERIFICATION_LOG.md
├── notes/
│   ├── failures.md              # What went wrong and why
│   ├── decisions.md             # Design decisions and rationale
│   └── approximations.md        # Master list of all approximations
├── output/                      # Rendered images
└── Project.toml                 # Julia project
```

---

## Phase 1: Literature Acquisition and Semantic Graph

**This phase MUST be completed before writing any simulation code.**

### 1.1 Download all required references

Obtain PDFs for each of the following. Store in `refs/` with systematic filenames. Record in `refs/INDEX.md`.

**Core physics — MUST HAVE:**

| Topic | Primary reference | What we need from it |
|-------|------------------|---------------------|
| Magnetised plasma EOS (fully ionised) | Potekhin & Chabrier, A&A 374, 213 (2001) | Eqs for free energy, pressure, entropy of e-p plasma in strong B |
| Magnetised plasma EOS (partial ionisation) | Potekhin & Chabrier, A&A 428, 787 (2004) | H ionisation equilibrium in strong B, occupation probabilities |
| Free-free opacity in strong B | Potekhin & Chabrier, A&A 399, 1007 (2003) | Appendix A: explicit ff absorption cross-sections, both modes |
| Magnetic Thomson/Compton scattering | Mészáros, "High-Energy Radiation from Magnetized Neutron Stars" (1992) ISBN 0-226-52094-4 | Ch 3–5: polarisation-dependent cross-sections, normal modes |
| Magnetic opacity — modern treatment | Potekhin, Chabrier & Ho, A&A 572, A69 (2014) | Updated opacities including vacuum polarisation, partial ionisation |
| Non-magnetic NS atmosphere | Heinke et al., ApJ 644, 1090 (2006) [NSATMOS] | Reference spectra for B=0 verification |
| Magnetic NS atmosphere models | Suleimanov, Potekhin & Werner, A&A 500, 891 (2009) | Full numerical scheme, published spectra for verification |
| Magnetic atmosphere review | Potekhin, Ho & Chabrier, JPCS 2016 / arXiv:1605.01281 | Review of current status, figure comparisons |
| Plane-parallel RT numerical methods | Mihalas, "Stellar Atmospheres" (2nd ed, 1978) ISBN 0-7167-0359-9 | Ch 6: Feautrier method, ALI, temperature correction |
| NS atmosphere RT (non-magnetic, open source) | Haakonsen et al., ApJ 749, 52 (2012) [McPHAC paper] | Verification benchmark, numerical scheme description |
| Kompaneets equation / Compton scattering | Rybicki & Lightman, "Radiative Processes in Astrophysics" ISBN 0-471-82759-2 | Ch 7: Compton scattering, Kompaneets equation |
| GR ray tracing — Schwarzschild | Pechenick, Ftaclas & Cohen, ApJ 274, 846 (1983) | Exact expressions for image of NS in Schwarzschild, elliptic integrals |
| GR ray tracing — modern treatment | Beloborodov, ApJ 566, L85 (2002) | Simplified approximate + exact formulae, benchmarks |
| GR ray tracing — Kerr | Dexter & Agol, ApJ 696, 1616 (2009) | Kerr geodesic integration |
| NS mass-radius / TOV equation | Shapiro & Teukolsky, "Black Holes, White Dwarfs, and Neutron Stars" ISBN 0-471-87316-0 | Ch 5–6: TOV equation, NS structure |
| TOV + EOS tables | Potekhin et al., A&A 560, A48 (2013) [BSk EOS] | Analytical EOS fits, Fortran routines described |
| Surface temperature from crustal transport | Greenstein & Hartke, ApJ 271, 283 (1983) | T(θ) for dipole B-field, analytic model |
| Surface temperature — modern | Potekhin, Pons & Page, SSRv 191, 239 (2015) / arXiv:1507.06186 | Review of thermal evolution, envelope models |
| Hydrogen in strong B — atomic physics | Ruder et al., "Atoms in Strong Magnetic Fields" (1994) ISBN 3-540-57499-0 | Energy levels, oscillator strengths, wavefunctions |
| Vacuum polarisation / QED birefringence | Ho & Lai, MNRAS 327, 1081 (2001) | Vacuum resonance, mode conversion |
| CIE 1931 colour matching functions | CIE 015:2004 or Stockman & Sharpe (2000) | Tabulated x̄(λ), ȳ(λ), z̄(λ) |
| XYZ → sRGB conversion | IEC 61966-2-1:1999 | Exact matrix, gamma curve |
| Tone mapping | Reinhard et al., SIGGRAPH 2002 | Photographic tone reproduction |
| NS cooling code reference | Page, "NSCool" documentation | Cross-check thermal evolution |
| X-PSI ray tracing reference | Riley et al., JOSS 8(82), 4977 (2023) | Verification of GR ray tracing implementation |
| McPHAC source code | GitHub: McPHAC/McPHAC | Structural reference + B=0 verification target |
| Compact object textbook | Shapiro & Teukolsky (1983) | NS interior, surface physics, GR effects |
| Radiative transfer in strong B | Lai, Rev. Mod. Phys. 73, 629 (2001) | Matter in strong B-fields review |

**Download strategy:**
- Use `arxiv.org/pdf/XXXX.XXXXX` for arXiv papers
- Use `doi.org/10.XXXX` and follow to publisher for journal versions
- For books: check if author has made chapters available; otherwise note ISBN and page ranges needed
- Clone McPHAC from GitHub: `git clone https://github.com/McPHAC/McPHAC.git refs/code/McPHAC`
- For CIE data: download from `cvrl.ioo.ucl.ac.uk` (Colour & Vision Research Laboratory — the canonical source for colour matching functions)

After downloading, run `sha256sum` on every file and record in `refs/INDEX.md`.

### 1.2 Build the semantic graph with `af`

Use the `af` command-line tool to construct a semantic graph of the entire physics pipeline. The graph has three types of nodes:

**CLAIM nodes:** A physics statement that can be true or false.
```
af claim add "The emergent specific intensity from a magnetised NS atmosphere depends on photon frequency, emission angle, local temperature, local B-field strength, and polarisation mode"
  --source potekhin_chabrier_2003
  --equation "n/a (conceptual)"
  --status accepted
```

**EQUATION nodes:** A specific mathematical relation tied to a ground-truth source.
```
af equation add "magnetic_ff_opacity_xmode"
  --latex "\kappa_{\rm ff}^{(X)} = ..."
  --source potekhin_chabrier_2003
  --equation_number "A.3"
  --page 1015
  --verified false
```

**APPROXIMATION nodes:** A simplification with named scope of validity.
```
af approx add "ground_landau_level_only"
  --description "Assume all electrons occupy the ground Landau level. Valid when ℏω_c >> kT, i.e., B >> 4.7e8 * (T/1e6)^2 G"
  --source meszaros_1992
  --page 47
  --affects [magnetic_ff_opacity_xmode, magnetic_ff_opacity_omode, magnetic_scattering]
  --validity_range "B > 1e10 G and T < 1e7 K"
```

**Required graph structure:** Build nodes and edges for at minimum:

1. **EOS subgraph:** Ideal gas → Coulomb corrections → Landau quantisation → partial ionisation (if included). Each correction is an APPROXIMATION node. Each formula is an EQUATION node.

2. **Opacity subgraph:**
   - Free-free (X-mode): equation, source, limits tested
   - Free-free (O-mode): equation, source, limits tested
   - Thomson scattering (X-mode): equation, source
   - Thomson scattering (O-mode): equation, source
   - Compton correction: equation, source
   - Bound-free (if partial ionisation): equation, source
   - Vacuum polarisation mode conversion (if B > 10^13 G): equation, source

3. **Radiative transfer subgraph:**
   - Feautrier discretisation: equation, source (Mihalas)
   - Boundary conditions (surface: no incoming; depth: diffusion approx): equations, sources
   - Temperature correction method: equation, source
   - Convergence criterion: definition, source

4. **GR ray tracing subgraph:**
   - Schwarzschild geodesic equation: equation, source
   - Impact parameter → surface coordinates map: equation, source
   - Redshift factor: equation, source
   - Solid angle element (for flux calculation): equation, source
   - Lensing of visible surface fraction: equation, source

5. **Surface model subgraph:**
   - Dipole B-field geometry: equation, source
   - T(θ) from Greenstein-Hartke or equivalent: equation, source
   - TOV for (M, R) from EOS: equations, source

6. **Colorimetry subgraph:**
   - Spectral integration with CIE functions: equation, source
   - XYZ → linear sRGB matrix: exact values, source
   - sRGB gamma curve: equation, source
   - Tone mapping operator: equation, source

Every edge in the graph must be typed: DEPENDS_ON, APPROXIMATES, VERIFIED_BY, IMPLEMENTS, DERIVED_FROM.

### 1.3 Literature review output

Before proceeding to Phase 2, produce:
- `refs/INDEX.md`: Complete index of all downloaded references with checksums
- `notes/approximations.md`: Master list of ALL approximations, each with name, validity range, source, and impact assessment (how much would the result change if this approximation were lifted?)
- `notes/decisions.md`: Record of all design decisions (e.g., "Start with fully ionised H atmosphere; partial ionisation deferred to Phase 4")
- Graph visualisation: `af` export of the full semantic graph

---

## Phase 2: Tracer Bullet — End-to-End with Simplified Physics

**Goal:** Get photons from the NS surface onto an image plane, with real (if simplified) physics at every step. A working pipeline with known approximations beats a perfect opacity module with no pipeline.

### 2.1 TOV Solver → (M, R, g)

Input: Central density ρ_c, EOS table (use BSk24 analytical fits from Potekhin et al. 2013).
Output: M, R, surface gravity g = GM/(R²√(1-2GM/Rc²)).

**Ground truth check:** Reproduce Table 1 of Potekhin et al. (2013) — M and R for given EOS to 4 significant figures.

### 2.2 Schwarzschild Ray Tracer

Input: Observer distance D (can be set to ∞ for parallel rays), image plane resolution, (M, R).
Output: For each pixel — surface coordinates (θ, φ) of hit point, emission angle cos(θ_e), gravitational redshift (1+z), or flag "misses star".

Implement using the elliptic integral formulae from Pechenick et al. (1983) or Beloborodov (2002). Start with Schwarzschild (no rotation).

**Ground truth checks:**
- Visible fraction of surface → must approach (1/2)(1 + √(1 - 2GM/Rc²))⁻¹ ... actually the exact formula: the entire surface with θ_max = π - arcsin(b_crit/R × √(1-2GM/Rc²)) where b_crit is from Pechenick et al.
- For M → 0 (flat space): visible fraction → exactly 1/2
- Reproduce Fig 1 of Pechenick et al. (1983)
- Compare against X-PSI's ray tracing if available

### 2.3 Surface Model

Input: Magnetic dipole parameters (B_pole, obliquity α), T_eff prescription.
Output: For each surface point — T_local, B_local, angle between B and surface normal.

Start with Greenstein-Hartke analytic model:
T(θ_B) = T_eq + (T_pole - T_eq) cos²(θ_B)
where θ_B is the magnetic colatitude.

Dipole field: B_r = B_pole cos(θ_B), B_θ = (B_pole/2) sin(θ_B), |B| = (B_pole/2)√(1 + 3cos²θ_B).

**Ground truth check:** Against Greenstein & Hartke (1983) equations.

### 2.4 Atmosphere: Blackbody Placeholder (then upgrade)

For the tracer bullet, start with a **modified blackbody**:
I_ν(ν, T, θ_e) = f_col⁴ × B_ν(ν, f_col × T) × (cos θ_e)^p

where f_col is a colour correction factor (~1.5–2.0, from Suleimanov et al.) and p is a limb-darkening exponent. This is a well-known approximation (hardness ratio fits) used in observational NS work.

**Record this as an APPROXIMATION node** with a clear flag: "PLACEHOLDER — to be replaced by self-consistent RT solver in Phase 3."

**Ground truth check:** In the limit f_col → 1, p → 0, must recover isotropic blackbody.

### 2.5 Redshift + Doppler

Apply relativistic transformations to each ray:
I_ν^obs = (1+z)^{-3} × I_ν'(ν' = (1+z)ν, θ_e)

If rotation: include Doppler boost δ = 1/[γ(1 - β·n̂)] where β is the surface velocity.

**Ground truth check:** Total observed luminosity must equal (1+z)^{-2} × (emitted luminosity seen by local observer). This is a consequence of I_ν/ν³ being a relativistic invariant.

### 2.6 Spectral → Colour

At each pixel, you now have I_ν^obs(ν) — the observed specific intensity as a function of frequency.

1. Integrate against CIE 1931 2° observer colour matching functions:
   X = ∫ I_ν^obs(ν) × x̄(ν) dν  (and similarly Y, Z)

2. XYZ → linear sRGB via the exact IEC 61966 matrix:
   [R_lin, G_lin, B_lin] = M × [X, Y, Z]

3. Apply sRGB gamma:
   C_srgb = 12.92 × C_lin  if C_lin ≤ 0.0031308
   C_srgb = 1.055 × C_lin^(1/2.4) - 0.055  otherwise

4. Tone map (the dynamic range is enormous — poles vs equator can differ by orders of magnitude in optical flux):
   Use Reinhard global operator: L_d = L / (1 + L) with appropriate key value.
   Or: offer multiple tone mapping choices (linear, Reinhard, filmic, ACES).

**Ground truth checks:**
- A 5778K blackbody (the Sun) must render as approximately (255, 249, 231) in sRGB after appropriate exposure.
- Verify CIE functions reproduce published chromaticity coordinates for standard illuminants (D65 etc).
- XYZ → sRGB matrix must match IEC 61966-2-1 exactly (to machine precision).

### 2.7 Tracer Bullet Integration Test

Render a 512×512 image of a canonical NS (M=1.4 M_☉, R=12 km, B_pole=10^12 G, T_pole=10^6 K, T_eq=5×10^5 K, distance 100 pc) viewed at inclination 60° to the magnetic axis.

Expected result: A dim blue-white disk, slightly larger than the geometric angular size (due to lensing), with faintly brighter polar caps. In false colour (mapping X-ray to visible), the structure should be dramatically visible.

**Output both:**
- True human-vision rendering (will be underwhelming — that's physically correct)
- False-colour multi-band rendering (map 0.1–0.5 keV → red, 0.5–2 keV → green, 2–10 keV → blue, or similar)

Save per-pixel spectral data (not just RGB) for later re-rendering with improved atmosphere models.

---

## Phase 3: Real Atmosphere — Self-Consistent Radiative Transfer

Replace the blackbody placeholder with a proper atmosphere calculation.

### 3.1 Non-Magnetic Atmosphere (B = 0)

Implement a Feautrier solver for a plane-parallel pure-H fully-ionised atmosphere in radiative equilibrium. This is the McPHAC problem. Use McPHAC's source code as a structural reference.

**Physics ingredients:**
- Free-free opacity: Kramers formula with Gaunt factor corrections. Source: Rybicki & Lightman Ch 5.
- Thomson scattering (anisotropic): Source: Chandrasekhar (1960) or McPHAC's treatment.
- Compton correction: Kompaneets equation. Source: Rybicki & Lightman Ch 7.
- Hydrostatic equilibrium: dP/dy = g (trivial in plane-parallel).
- Radiative equilibrium: ∫_0^∞ (κ_ν J_ν - η_ν) dν = 0 at each depth.

**Ground truth checks:**
- Reproduce McPHAC spectra to <1% across 0.01–10 keV for at least 3 different (T_eff, log g) combinations.
- Verify Eddington-Barbier approximation in Rayleigh-Jeans limit.
- Verify limb darkening law against Chandrasekhar's H-function in the grey atmosphere limit.

### 3.2 Magnetic Atmosphere

Extend to two polarisation modes in a magnetised plasma. This is the central physics challenge.

**New ingredients:**
- Magnetic free-free opacities: Potekhin & Chabrier (2003) Appendix A. Both X-mode and O-mode.
- Magnetic scattering cross-sections: Mészáros (1992) Ch 5.
- Two coupled RT equations (or uncoupled if ignoring mode conversion).
- Modified hydrostatic equilibrium (magnetic pressure gradient is negligible for B < 10^15 G — verify this claim and record as APPROXIMATION).

**Implementation order:**
1. Implement magnetic opacities as standalone functions with extensive unit tests.
2. Test against B→0 limit: must recover non-magnetic opacities.
3. Test against known values from Potekhin & Chabrier (2003) tables (these tables are derived from the published equations — reproduce the equations, then check your code against the tables as a verification step, not as a dependency).
4. Plug into the Feautrier solver (now with two modes).
5. Run to convergence for a canonical case (T_eff = 10^6 K, log g = 14.3, B = 10^12 G).
6. Compare emergent spectrum against Suleimanov et al. (2009) published figures.

### 3.3 Atmosphere Table Generation

Once the solver works, pre-compute a grid of atmosphere models:
- T_eff: 10^{5.5} to 10^{7.0} in steps of 0.1 dex (16 points)
- log g: 13.5 to 14.8 in steps of 0.2 (7 points)
- B: 0, 10^{10}, 10^{11}, 10^{12}, 10^{13}, 10^{14} G (6 points)
- cos θ_e: 0.05 to 1.0 in 20 points

Each model outputs I_ν(E, cos θ_e) for both polarisation modes, on an energy grid from 0.01 eV to 100 keV (covering optical through hard X-ray).

**These are YOUR tables, generated by YOUR code, from published equations.** They are fully reproducible. Store them in a documented format (HDF5 or FITS with complete metadata including code version, git hash, and all input parameters).

---

## Phase 4: Extensions (after tracer bullet works)

In priority order:

### 4.1 Partial Ionisation
Add bound states using Potekhin & Chabrier (2004) formalism. This matters for T_eff < 5×10^5 K and B > 10^12 G. Adds bound-free opacity and modifies EOS.

### 4.2 Vacuum Polarisation and Mode Conversion
Implement the vacuum resonance where the plasma and vacuum contributions to the dielectric tensor cancel. Photons passing through this layer can convert between X and O modes. Source: Ho & Lai (2001). Matters for B > 10^13 G.

### 4.3 Condensed Surface
At very high B and/or low T, the surface may be a condensed solid/liquid rather than a gaseous atmosphere. Emission from a condensed surface: van Adelsberg et al. (2005). This is an alternative to the atmosphere model, not an addition.

### 4.4 Kerr Ray Tracing
Add rotation. Replace Schwarzschild geodesics with Kerr. Source: Dexter & Agol (2009). Adds Doppler boosting and frame dragging.

### 4.5 Magnetosphere
Volume rendering of pair plasma in the magnetosphere. This is where Lyr.jl's MC volume rendering engine becomes essential. The magnetosphere is optically thin in most cases but contributes resonant cyclotron scattering features.

### 4.6 Interstellar Medium Absorption
For "realistic" appearance: apply interstellar absorption using Wilms, Allen & McCray (2000) cross-sections. This reddens the spectrum and creates absorption edges. Only relevant if visualising a specific known NS at a known distance with known N_H.

---

## Technical Requirements

### Julia Packages (install as needed)
```julia
# Core numerics
DifferentialEquations    # TOV solver, possibly RT integration
QuadGK                   # Numerical integration (spectral integrals)
Interpolations           # Atmosphere table interpolation
SpecialFunctions         # Bessel, Laguerre, elliptic integrals
Elliptic                 # Complete and incomplete elliptic integrals (for geodesics)
StaticArrays             # Performance-critical small arrays
LinearAlgebra            # LAPACK for tridiagonal solves (Feautrier)

# Data I/O
HDF5 / FITS / JLD2       # Atmosphere table storage
FileIO                   # Image output

# Visualisation
Colors                   # Colour space conversions
Images                   # Image construction
Makie / CairoMakie       # Plotting and visualisation

# Verification
Test                     # Unit testing
BenchmarkTools           # Performance
```

### Non-Julia Dependencies (acceptable for Phase 2 tracer bullet)
- McPHAC (C code): Clone for verification. Do NOT depend on it at runtime.
- CFITSIO (C library): Only if reading FITS files for verification against XSPEC tables.
- Potekhin Fortran routines: Download for verification. Wrap with `ccall` if needed for cross-checking, but the Julia implementation must be independent.

### Coding Standards
- Every physics function gets a docstring citing its source: `@doc "Free-free opacity, X-mode. Potekhin & Chabrier (2003) Eq. A.3" function ...`
- CGS units throughout the physics code (this is astrophysics, fight me). Define a `PhysicalConstants` module with CODATA 2018 values, each citing the source.
- No magic numbers. Every constant gets a name and a source citation.
- ~200 lines per file maximum. Decompose aggressively.
- Test files mirror source files: `src/opacity/magnetic_ff.jl` → `test/opacity/magnetic_ff_test.jl`

---

## Verification Matrix

Every module must pass ALL checks before being marked complete.

| Module | Limit Check | Literature Check | Cross-Code Check |
|--------|-------------|-----------------|-----------------|
| TOV solver | M→0: Newtonian limit | Potekhin+(2013) Table 1 | Against NSCool if possible |
| EOS (B=0) | Ideal gas limit | Potekhin+(2001) published values | — |
| EOS (B>0) | B→0: recover non-mag EOS | Potekhin & Chabrier (2001) table values | Against Potekhin Fortran |
| Opacity ff (B=0) | Kramers opacity at low ν | Rybicki & Lightman examples | Against McPHAC |
| Opacity ff (B>0) | B→0: recover non-mag opacity | Potekhin & Chabrier (2003) Fig 1 | Against Potekhin Fortran |
| Opacity scat (B=0) | Thomson cross-section | σ_T = 6.652e-25 cm² | Against McPHAC |
| Opacity scat (B>0) | ν≪ν_c: suppressed X-mode | Mészáros (1992) Ch 5 examples | — |
| RT solver (B=0) | Grey atmosphere: Hopf function | McPHAC spectra (<1%) | McPHAC direct comparison |
| RT solver (B>0) | B→0: recover non-mag spectra | Suleimanov+(2009) Figs 4-6 | Against XSPEC nsmax tables |
| Ray tracer | M→0: flat space projection | Pechenick+(1983) Fig 1 | Against X-PSI if possible |
| Colorimetry | 5778K BB → solar white | CIE daylight locus | — |
| Tone mapping | Identity for low dynamic range | Reinhard+(2002) Fig 7 | — |

---

## Exit Criteria

The project is "done" (for a given phase) when:

**Phase 1:** All references downloaded, checksummed, indexed. Semantic graph has ≥50 equation nodes, ≥20 approximation nodes, all linked to local PDFs. `notes/approximations.md` is complete.

**Phase 2:** Tracer bullet renders a 512×512 image. All verification checks in the matrix above pass for the simplified (blackbody) physics. Per-pixel spectral data is saved. Both true-colour and false-colour renders are produced.

**Phase 3:** Self-consistent non-magnetic atmosphere reproduces McPHAC to <1%. Magnetic atmosphere reproduces Suleimanov et al. (2009) to <5% (accounting for differences in input physics). Full re-render with proper atmosphere.

**Phase 4+:** Extensions as prioritised above.

---

## Anti-Patterns — Things That Must Never Happen

1. **Implementing an equation "from memory"** without finding it in a local PDF. Even if you "know" the Planck function, cite Rybicki & Lightman Eq 1.51 (or wherever) and record it.

2. **Downloading a table from a URL and using it as an input** without being able to regenerate it from published equations. Tables are for VERIFICATION, not for PRODUCTION.

3. **Declaring a test "passed" without printing the comparison values.** Always print: expected value (with source), computed value, relative error, tolerance.

4. **Silent fallbacks.** If an interpolation extrapolates, if a solver doesn't converge, if a value is out of range — it MUST be logged. Never silently clamp, extrapolate, or substitute a default.

5. **"It looks right."** Visual inspection is not a verification strategy. Every claim of correctness must be backed by a quantitative comparison against a published value.

6. **Skipping the semantic graph.** The graph is not bureaucracy. It is how we know that every piece of the pipeline is grounded. If you can't trace a number back to a local PDF through the graph, the number is not trustworthy.

---

## Getting Started — First Session Checklist

1. Create the directory structure above.
2. Begin downloading papers. Start with the arXiv papers (fastest). For each:
   - `wget https://arxiv.org/pdf/XXXX.XXXXX -O refs/author_year_journal.pdf`
   - `sha256sum refs/author_year_journal.pdf >> refs/CHECKSUMS`
   - Add entry to `refs/INDEX.md`
3. Clone McPHAC: `git clone https://github.com/McPHAC/McPHAC.git refs/code/McPHAC`
4. Initialise Julia project: `julia -e 'using Pkg; Pkg.generate("NeutronStar")'`
5. Begin `af` graph construction with the highest-level claims and work downward.
6. Implement `src/constants.jl` — all physical constants with CODATA citations.
7. Implement TOV solver (simplest self-contained module with clear verification target).
8. **Do not** start on the atmosphere code until Phase 1 is substantially complete.

---

*This prompt is a living document. Update it as decisions are made and recorded in `notes/decisions.md`.*
