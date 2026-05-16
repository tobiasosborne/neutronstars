# Architecture Review — NeutronStar.jl

Reviewer: architecture lens only (correctness vs paper equations, test coverage, tidiness handled elsewhere).
Date: 2026-05-16, commit `95c76b3`.

## Executive summary

- **The renderer silently uses the wrong physics.** `AtmosphereGrid` is built with a hardcoded `θ_B = 0.0` (pole-on field, see `src/pipeline/atmosphere_grid.jl:79`) but the renderer asks for spectra at arbitrary surface points without re-passing `θ_B`. Every "magnetic atmosphere" render in the project is in fact a pole-only atmosphere mapped onto an oblate magnetic geometry. This is undocumented, the type system does not catch it, and `lookup_spectrum` is silent about it.
- **The exported public API contains a latent EOS bug.** `AtmosphereStructure.update_atmosphere!` (line 144) drops the factor of 2 from the fully-ionised-H mean molecular weight (`ρ = m_p P / (k_B T)` instead of `/(2 k_B T)`). The function is exported, documented, and dead — the iteration loop uses a private copy (`RTAtmosphere._update_structure!`) that has the right formula. Any user who follows the public API will get a 2× density error.
- **No provenance, no parallelism, no configuration layer.** `AtmosphereSpectrumGrid` is just a struct of arrays with no record of how it was built (K, M, N, tolerances, gaunt-table hash, code version). The `_compute_magnetic_opacities!` and per-frequency Feautrier loops are pure serial despite being embarrassingly parallel. There is no logging abstraction (39 `println`/`@printf` sites). There is no `Configuration` object — solver knobs (`tol`, `flux_damping`, `y_max`, `N_eff` cutoffs, 80, 10, 0.265) are scattered as defaults and magic numbers across five files.

---

## Findings

### 1. Renderer silently drops magnetic-field obliquity (`θ_B`)
- **Where:** `src/pipeline/atmosphere_grid.jl:79`; `src/pipeline/render.jl:300-313`; `src/pipeline/atmosphere_grid.jl:94-137` (`lookup_spectrum`).
- **What:** The grid axes are `(T_eff, B)` only; `build_atmosphere_grid` passes a constant `θ_B = 0.0` to `solve_magnetic_atmosphere`. The renderer computes `θ_B = magnetic_colatitude(...)` per pixel but throws it away — `lookup_spectrum(grid, T_local, B_local, ν_grid_em, cos_θe)` has no `θ_B` parameter. The grid struct therefore represents one specific angle (B parallel to surface normal). For a dipole geometry, `θ_B` varies from 0 at the magnetic pole to π/2 at the magnetic equator and is the dominant variable for X/O mode opacity ratio (Suleimanov+ 2009 Fig 2).
- **Why it matters:** Every "magnetic" rendering in `output/` (RX J1856 rotation, spectral sweep) silently uses pole-on opacities everywhere. The hot-spot anisotropy people see is purely from the Greenstein–Hartke T(θ) map, not from RT. Validating the renderer against observed pulsed-fraction shapes will give a false-positive agreement and then a much worse mismatch when `θ_B` is wired in. The HANDOFF.md "Working And Verified" claim for Phase 3c is materially misleading.
- **Severity:** critical.
- **Suggested fix:** Add a `θ_B_grid::Vector{Float64}` axis to `AtmosphereSpectrumGrid`; build with at least `θ_B ∈ {0, π/4, π/2}` initially; route `θ_B` through `lookup_spectrum` and `render_spectral_cube`. Until that lands, rename current artifact `AtmosphereSpectrumGridPoleOn` so future agents cannot mistake it for a full table.

### 2. Public `update_atmosphere!` has the wrong EOS factor; the iteration loop bypasses it
- **Where:** `src/atmosphere/atm_structure.jl:138-164` (exported, documented). `src/atmosphere/rt_atmosphere.jl:116, 162-184` calls a private `_update_structure!` instead.
- **What:** `update_atmosphere!` line 144: `col.ρ[i] = m_p * col.P[i] / (k_B * col.T[i])` — missing factor of 2. The hydrostatic build (`build_atmosphere`, line 93) uses `/ (2 k_B T)`; `_update_structure!` line 165 also uses `/ (2 k_B T)`. The bad function is imported by `RTAtmosphere` (line 21) and re-exported through `NeutronStar.jl` line 77 but never actually called.
- **Why it matters:** A user who writes their own iteration loop using only the exported API will get densities 2× too high, opacities 2× too high, optical depth 2× too high, and a converged (but wrong) atmosphere. Dead public code is worse than no public code. Also indicates that the AtmosphereColumn invariant "ρ corresponds to (P,T) under ionised-H EOS" is not enforced anywhere.
- **Severity:** high.
- **Suggested fix:** Delete `update_atmosphere!` (or fix the formula AND make `RTAtmosphere._update_structure!` call it, removing the duplicate). Better: encapsulate the EOS as `AtmosphereStructure.density_from_PT(P, T)` and call it from both sites and from `_build_initial_column`.

### 3. `AtmosphereSpectrumGrid` has no provenance metadata
- **Where:** `src/pipeline/atmosphere_grid.jl:28-37`.
- **What:** Fields: `T_grid, B_grid, g_s, ν_grid, μ_grid, I_cache`. No record of: K, M, N, max_iter, tol, anisotropic flag, Gaunt-table path/hash, code version, build date, mode coupling assumptions, θ_B value, `flux_damping`, or convergence status per grid cell. `SpectralImageCube.atmosphere_id` (`src/pipeline/render.jl:245`) is a hardcoded literal `"AtmosphereGrid_v1"`.
- **Why it matters:** When the HDF5 caching listed in HANDOFF.md backlog ships, you'll save grids and then have no way to know whether a stored grid corresponds to the current code. Two atmospheres built with different `flux_tol` will look identical on disk. Reproducibility is gone, and reviewers cannot answer "is this rendered image current?".
- **Severity:** high.
- **Suggested fix:** Add `struct AtmosphereGridProvenance { code_sha::String, build_time::DateTime, K, M, N, max_iter, tol_T, tol_flux, flux_damping, gaunt_sha::String, theta_B::Float64, converged::Matrix{Bool}, iters::Matrix{Int} }` and embed it in the grid struct and in any saved file.

### 4. No `θ_B` parameter in solver-result struct, even though the angle drives the physics
- **Where:** `src/atmosphere/magnetic_atmosphere.jl:32-45` (`MagneticAtmosphereResult`) records `B` and `θ_B`. Good. But `AtmosphereSpectrumGrid` does not — the angle disappears at the grid layer. Combined with finding #1, the grid completely loses the angle dimension.
- **Why it matters:** Anybody trying to plot pulse profiles for a non-zero magnetic obliquity will have no way to do it without rebuilding the entire grid abstraction.
- **Severity:** high (subordinate to #1).
- **Suggested fix:** See #1.

### 5. `lookup_spectrum` interpolates poorly: nearest-neighbour in ν, linear in μ, linear in T/B
- **Where:** `src/pipeline/atmosphere_grid.jl:101-137` and `144-163`.
- **What:** Line 127 `kν = _nearest(grid.ν_grid, ν_obs[iν])` — nearest-neighbour, no log-log interpolation. `_nearest` (line 194) uses `argmin(abs.(grid .- x))` which **allocates** a fresh array of size `length(grid)` for every call. A 512×512×50 cube does this 13M times. The dedicated `rt_emergent_spectrum` (rt_atmosphere.jl:191) does log-log frequency interpolation and bilinear angle interpolation — strictly better, and not used by the grid.
- **Why it matters:** (a) stair-stepped spectra around emission lines and cyclotron features, (b) the renderer is far slower than necessary, (c) two pieces of nearly-identical interpolation logic must be kept in sync. The unused `lookup_spectrum_scalar` is a third copy.
- **Severity:** medium.
- **Suggested fix:** Replace `_nearest` with `searchsortedlast` + log-log interpolation; drop `lookup_spectrum_scalar` or define it as a wrapper around the vector version; consider extracting interpolation to a small module both the grid and `rt_emergent_spectrum` can call.

### 6. Magnetic-emergent-spectrum interpolation is duplicated in a script
- **Where:** `scripts/render_rxj1856_visible_magnetic_atmosphere.jl:24-62`. Compare `src/atmosphere/rt_atmosphere.jl:191-242` (`rt_emergent_spectrum`).
- **What:** `rt_emergent_spectrum` only handles the non-magnetic 2-D `I_emergent[k, j]`. For magnetic results (`K × M × 2`), the script reimplements the same logic with an inner `sum(... for mode in 1:2)`. There is no exported `magnetic_emergent_spectrum`.
- **Why it matters:** Any future analysis or downstream code that wants to read magnetic emergent spectra will reinvent this. The fact that one of the two existing copies lives in `scripts/` (not under `src/`) means it can't be tested.
- **Severity:** medium.
- **Suggested fix:** Promote the script's logic to `MagneticAtmosphere.magnetic_emergent_spectrum` and export it. Even better: make `I_emergent` a typed wrapper that knows whether it has a polarization axis, and dispatch.

### 7. `MagneticAtmosphereResult`'s `I_emergent` axis order is `K × M × 2`; the non-magnetic version is `K × M`. The renderer flattens them through ad-hoc summation
- **Where:** `src/atmosphere/magnetic_atmosphere.jl:42, 171`; `src/atmosphere/rt_atmosphere.jl:128-131`; `src/pipeline/atmosphere_grid.jl:74, 82`.
- **What:** Non-magnetic `I_emergent[k, j]` (K×M). Magnetic `I_emergent[k, j, mode]` (K×M×2). `build_atmosphere_grid` collapses the magnetic case via `r.I_emergent[:, :, 1] .+ r.I_emergent[:, :, 2]` into a K×M matrix and stores it in `I_cache::Array{Matrix{Float64}, 2}` — so the grid loses all polarization information.
- **Why it matters:** (a) Polarized renders (HANDOFF.md mentions phase-4 IXPE-like products) are impossible from the current cache. (b) The collapsed format means observers integrating over modes implicitly assume both modes escape with the same redshift / angle dependence at the renderer, which is true in Schwarzschild but not in Kerr where the polarization basis rotates along the geodesic. (c) The shape contract is enforced by neither type signature; a user could pass a 4-D result and the broadcast would silently mis-sum.
- **Severity:** medium (now), high (when polarization or Kerr lands).
- **Suggested fix:** Introduce `struct EmergentIntensity{N}` where `N ∈ {2, 3}` is the polarization axis count; have `lookup_spectrum` return a polarized spectrum; provide an explicit `sum_modes` adapter.

### 8. Massive coupling between rt_atmosphere ↔ atm_structure: temp-correction & opacities reach into `AtmosphereColumn` fields directly
- **Where:** `src/atmosphere/rt_atmosphere.jl:162-184`, `src/atmosphere/temp_correction.jl:50-228`, `src/atmosphere/feautrier.jl:50-101`.
- **What:** `_update_structure!`, `_build_rybicki_system!`, `_solve_single_frequency` all index `col.T[i]`, `col.ρ[i]`, `col.P[i]`, `col.κ[i,k]`, `col.k_total[i,k]`, `col.ρ_alb[i,k]`, `col.τ[i,k]`, `col.y[i]`, `col.ν_grid[k]`, `col.k_R[i]`, `col.gaunt`, `col.N`, `col.K`, `col.σ_scat`. Effectively every field is a public API. The magnetic atmosphere solver does not use `AtmosphereColumn` at all — it carries `(y, T, ρ, ν_grid, κ, k_total, ρ_alb, τ)` as bare arrays passed by position through ~7 internal functions.
- **Why it matters:** (a) The two RT pipelines (magnetic, non-magnetic) cannot share temp-correction or Feautrier code because they use incompatible data formats. (b) Refactoring `AtmosphereColumn` (e.g. adding `θ_B`) ripples through every interior function. (c) Two implementations of `_tridiag_lu_forward` / `_tridiag_lu_back` exist (`temp_correction.jl:234-269` and `magnetic_atmosphere.jl:607-627`) because of this.
- **Severity:** high.
- **Suggested fix:** Define a `RTState` interface (T, opacities, depths, frequencies, polarization axis, BCs) that both modules use. Move tridiagonal solver to a `Linalg` helper module. The view-based magnetic solver can then be a few lines on top of a generic Feautrier kernel.

### 9. Implicit array layout `[i, k, j]` vs `[k, j, mode]` is not documented anywhere checkable
- **Where:** `src/atmosphere/magnetic_atmosphere.jl:85-100` (`κ[N, K, 2]`, `P_all[N, M, K, 2]`, `I_emergent[K, M, 2]`); `src/atmosphere/feautrier.jl:58` (`P_all[N, M, K]`); `src/pipeline/atmosphere_grid.jl:34, 122` (`I_cache[iT, iB]` whose inner element is `K × M`).
- **What:** Four different array shape conventions. `P_all` is `[depth, angle, freq, mode]`, but `I_emergent` is `[freq, angle, mode]` (transposed first two). `I_cache` is a Matrix-of-Matrices where the outer dims are `(T, B)` (T first) and inner are `(K, M)`. The docstrings sometimes say "K × M" but the code is `zeros(N, M, K, ...)`. There is no type-level contract.
- **Why it matters:** Adding any new axis (e.g. `θ_B` or polarization mode) requires re-checking every nested loop in five files. Bug class: silent transpose. Currently the test suite has no shape-contract test.
- **Severity:** medium.
- **Suggested fix:** Wrap the multi-axis arrays in named tuples or small `AxisArray`-style structs (or at minimum add `@assert size(P_all) == (N, M, K)` at every boundary). Adopt one convention (e.g. always `(depth, angle, freq, [mode])`) and refactor `I_emergent` to match.

### 10. Hidden mutable state in `AtmosphereColumn`; iteration loops mutate fields in place but `MagneticAtmosphereResult` and `AtmosphereResult` only snapshot `copy(col.T)` and `copy(col.y)`
- **Where:** `src/atmosphere/atm_structure.jl:26` (`mutable struct AtmosphereColumn`); `src/atmosphere/rt_atmosphere.jl:137-138` (`copy(col.T), copy(col.y)`).
- **What:** All other fields (κ, ρ, P, τ, k_R) are not snapshotted, so a `AtmosphereResult` is meaningfully incomplete — you cannot reconstruct the converged opacity profile from it. A user holding an `AtmosphereResult` cannot re-run the Feautrier solve at a finer μ grid without re-iterating.
- **Why it matters:** Slot debugging and downstream analyses (line profiles, depth-of-formation) need the full structure. The non-result fields effectively die with the function call.
- **Severity:** medium.
- **Suggested fix:** Either snapshot the whole `AtmosphereColumn` into the result, or replace `mutable struct` with an immutable + functional updates. The latter also opens the door to parallel iteration over different T_eff at once.

### 11. No concurrency anywhere; embarrassingly-parallel loops are serial
- **Where:** `grep -rn 'Threads\.' src/` returns nothing.
- **What:** Hot loops where each iteration is independent:
  - `build_atmosphere_grid` over `(iT, iB)` (`src/pipeline/atmosphere_grid.jl:65-87`) — currently 12 minutes for the magnetic case per HANDOFF.md.
  - `_compute_magnetic_opacities!` over `(i, k)` (`src/atmosphere/magnetic_atmosphere.jl:276-298`) — both opacity calls are pure.
  - `_solve_feautrier_mode!` over `k` (`src/atmosphere/magnetic_atmosphere.jl:319-375`) — single-frequency solves are independent.
  - `solve_feautrier_all` over `k` (`src/atmosphere/feautrier.jl:63-98`) — same.
  - `render_spectral_cube` over `(i, j)` pixels (`src/pipeline/render.jl:292-319`).
- **Why it matters:** Modern desktops are 8–32 cores. Grid build is a 12× speedup hanging fruit. Render is similar. Phase 4 Kerr ray tracing will be 10–100× more expensive and will need thread or task parallelism baked in from the start.
- **Severity:** medium (today), high (Phase 4).
- **Suggested fix:** Add `Threads.@threads` to the outer grid build loop first (each cell is fully independent and allocates its own scratch). Audit any shared scratch buffers in `_solve_feautrier_mode!` before parallelising inside an iteration — currently `dtau`/`B_tilde`/`Q_tilde`/`C_store` are allocated per frequency, so it should be safe. Document thread-safety contracts on each solver entry point.

### 12. `argmin(abs.(grid .- x))` allocates a vector in the inner render loop
- **Where:** `src/pipeline/atmosphere_grid.jl:194-196`. Called inside `lookup_spectrum` per pixel per frequency.
- **What:** Allocates `length(grid)` Float64s each call. With K~50, 512² pixels, 50 cube frequencies, that's roughly 6.5 GB of garbage for a single render.
- **Why it matters:** GC pressure, cache pollution, blocking parallelism. Also masks the deeper issue of nearest-neighbour ν interpolation.
- **Severity:** medium.
- **Suggested fix:** Replace with `searchsortedlast` + log-log interp (also fixes finding #5).

### 13. The `Gaunt` table is mandatory in `solve_magnetic_atmosphere` even when B > 0 (where it is unused)
- **Where:** `src/atmosphere/magnetic_atmosphere.jl:56-65`; `_compute_magnetic_opacities!` at `src/atmosphere/magnetic_atmosphere.jl:272-298`.
- **What:** Signature is `solve_magnetic_atmosphere(T_eff, g_s, B, θ_B, gaunt::GauntTable; ...)`. The B>0 path uses `mode_absorption` / `mode_scattering` which have their own internal Coulomb logarithm and never touch the Gaunt table. The B=0 path uses `kappa_ff(ν, T, ρ, gaunt)`. The test at `test/runtests.jl:134` already exploits this loophole by passing `nothing` (which works because the type assertion is positional, not strict). The test would crash for B=0.
- **Why it matters:** API lies about its dependencies. A user wanting a pure magnetic run still has to load a 75kB Gaunt table from a hardcoded path (`refs/code/McPHAC/gffgu.dat`). The Gaunt-table file is also a build-time dependency outside the package directory.
- **Severity:** medium.
- **Suggested fix:** Make `gaunt::Union{Nothing, GauntTable}=nothing` and assert `gaunt !== nothing` only in the B=0 branch. Or factor the B=0 case into a separate function. Long-term: ship a small Gaunt table inside the package as an artifact and load it lazily.

### 14. Hardcoded reference paths in `Renderer`
- **Where:** `src/pipeline/render.jl:103` and `:344`: `cmfs_path = joinpath(dirname(dirname(@__DIR__)), "refs", "cvrl_cie1931_2deg.csv")`.
- **What:** The colour matching function path is computed from the source file location. If `NeutronStar.jl` is added as a registered package or vendored into another project, this path breaks. Also `3.0856e18` at line 325 duplicates the `pc` constant from `PhysicalConstants` (which is defined to higher precision as `3.0856775814913673e18`).
- **Why it matters:** Distribution friction; Phase 4 says "extract Lyr.jl as nested module" — package-relative file lookup is at exactly the wrong layer.
- **Severity:** low–medium.
- **Suggested fix:** Use Julia `Artifacts` for the CMF data; use `PhysicalConstants.pc` instead of literal.
- **Status (E13, 2026-05-16):** CMF half resolved — file vendored under `src/data/cvrl_cie1931_2deg.csv` and both load sites switched to `joinpath(@__DIR__, "..", "data", ...)`. Decided against `Artifacts.toml` for a 23 KB CSV: hosting + `LazyArtifacts` overhead is disproportionate, and `src/` content ships automatically with the package. The `3.0856e18` literal is still outstanding.

### 15. `_flux_ratio` diagnostic in `AtmosphereGrid` is numerically wrong
- **Where:** `src/pipeline/atmosphere_grid.jl:168-177`.
- **What:**
  ```julia
  F += 2π * μ[j] * I_em[k, j] * (μ[j] > 0 ? 2.0 : 0.0) * dν
  ```
  Missing the Gauss-Legendre weight `w[j]`. Since `gauss_legendre_half` returns only positive μ, the `μ[j] > 0` branch always picks 2.0 — meaning the function pretends to sum `2π·μ·I·2·dν` but should sum `2π·μ·I·w·dν` (or with an extra factor of 2 if `I=2P` substitution wasn't already done). The printed `F/σT⁴` is off by ~length(μ).
- **Why it matters:** This number is printed during grid construction and visible to anyone running the tutorial. It reassures the user that the grid converged when the underlying solver may have not (or vice versa). It also makes future agents distrust real flux numbers from `solve_atmosphere`.
- **Severity:** medium (silent diagnostic lying).
- **Suggested fix:** Either remove this and call the canonical `_bolometric_flux` from `RTAtmosphere`, or fix the formula and add a test that pins it.

### 16. Module include order is fragile and undocumented
- **Where:** `src/NeutronStar.jl:1-98`.
- **What:** The order is comment-driven: "blackbody must come before magnetic_modes" (line 35), "dielectric tensor must come before magnetic_modes" (line 31). There are 16 includes and a misordering will produce `UndefVarError`s during precompilation that point to the wrong file. There is no test that catches an incorrect include order short of full precompile.
- **Why it matters:** Adding any module (HDF5 cache, partial ionisation, Kerr tracer) requires manually picking the right insertion point. The dependency graph is implicit in the include comments rather than declared in code.
- **Severity:** low.
- **Suggested fix:** Either accept this is how Julia works and add a small ASCII module-graph diagram in a header comment, or move each module to its own file with `Pkg.develop`-style sub-packages so that `using` resolves the order automatically. The latter is overkill for the current scale.

### 17. `PotekhinTableReader` is defined but never included
- **Where:** `src/eos/potekhin_table_reader.jl` (235 lines, complete module). `grep` shows no `include` of it in `src/NeutronStar.jl` or anywhere else.
- **What:** Dead module. References say the Potekhin tables are used as a verification target, but the parser is orphaned.
- **Why it matters:** Either (a) the table reader should be wired into `verification/` and used to cross-check `MagneticModes.rosseland_magnetic` against Potekhin's tabulated `K_∥, K_⊥`, or (b) it should be deleted. Right now it's a maintenance liability — a refactor could break it silently and nobody would notice.
- **Severity:** low.
- **Suggested fix:** Either include + add a verification test that compares our Rosseland opacities against `refs/potekhin_tables/`, or delete.

### 18. No structured error types — every failure is `@assert` or `error("...")` with a string
- **Where:** ~40 `@assert` sites across `src/`. Examples: `src/eos/tov.jl:64-65` (`@assert ρ_c > 1e10`), `src/atmosphere/atm_structure.jl:74-75`, `src/atmosphere/magnetic_atmosphere.jl:65-66`. `src/eos/bsk_eos.jl:235` (`error("Bisection did not converge for P=$P after 200 iterations")`).
- **What:** Errors carry no machine-readable distinction between "user gave bad input" / "physics inconsistency" / "numerical convergence failure" / "table out of range". Callers (e.g. `mass_radius_curve` at `src/eos/tov.jl:215`) catch `Exception` generically.
- **Why it matters:** When you wire Kerr ray tracing on top of TOV in Phase 4, you'll want to know whether a failed solve is "central density too high" (skip this point in the M–R curve) vs "horizon reached" (probably a bug). Right now you can only string-match. Tests cannot assert `@test_throws CertainSpecificError`.
- **Severity:** low–medium.
- **Suggested fix:** Define `abstract type NSError <: Exception end` with `BadInputError`, `NumericalError`, `ConvergenceError`, `TableOutOfRangeError`. Cheap; pays off later.

### 19. No units handling — all CGS is convention only
- **Where:** Throughout. Docstrings list units; nothing enforces them.
- **What:** `NSParams` has `R_km::Float64` (km) but `T_pole::Float64` (K) — same type, different units. `g_s` is `cm/s²`. `B` is Gauss. A user passing `R = 12000` (meters? km? cm?) gets silent bad answers.
- **Why it matters:** This is a known long-term risk; Unitful.jl is the standard remedy. The Phase 4 magnetosphere volume renderer will introduce many more length scales and the room for a units bug grows.
- **Severity:** low (today), medium (when external collaborators arrive).
- **Suggested fix:** At minimum, define newtypes for the most error-prone quantities (`struct Density; ρ::Float64 end`). Or commit to Unitful and refactor the boundary functions only.

### 20. `verbose::Bool` is the only configuration channel; no structured logging
- **Where:** Every solver entry point. 39 `println`/`@printf` sites (grep above).
- **What:** Logging is on/off boolean; no log levels, no timestamps, no per-module filters. `@printf("  iter %2d: ...")` inside iteration loops cannot be redirected to a file without also losing the convergence messages.
- **Why it matters:** When the grid build takes 12 minutes, the user has no way to monitor progress without watching stdout live. Test runs print thousands of lines. Tooling like log aggregation/grafana is out of reach.
- **Severity:** low.
- **Suggested fix:** Replace `verbose &&` with `@debug`/`@info`/`@warn` from Julia's `Logging` stdlib; expose a single `with_logger(...)` wrapper for scripts.

### 21. Solver tolerances and physical thresholds are scattered magic numbers
- **Where:** `0.265 T_eff` (`atm_structure.jl:178`, `magnetic_atmosphere.jl:224`); `τ ≥ 80` (`feautrier.jl:67, 137`; `magnetic_atmosphere.jl:211, 343`); `τ_eff > 10` (`magnetic_atmosphere.jl:343`); `y_max = 1e5` (`magnetic_atmosphere.jl:190`); `flux_damping=0.5` (`magnetic_atmosphere.jl:63`); `0.3 / max_dT` damping cap (`rt_atmosphere.jl:93-97`, `magnetic_atmosphere.jl:131-134`); `[0.9, 1.1]` flux clamp (`magnetic_atmosphere.jl:126`).
- **What:** Each constant lives at its use site. The 0.265 boundary constant appears twice and would need to change in lockstep. The τ_max=80 appears four times. There is no single `SolverParameters` struct, no `MagneticSolverConfig`.
- **Why it matters:** Tuning these for Phase 4 partial-ionisation or higher T_eff regimes requires hunting through five files. Two of the four `τ ≥ 80` sites might silently fall out of sync.
- **Severity:** medium.
- **Suggested fix:** A `module SolverDefaults` exporting `const TAU_DIFFUSION = 80.0`, `const T_SURFACE_FRAC = 0.265`, etc., or a `SolverParams` struct passed through the call chain.

### 22. `n_e` definition uses `ρ / (m_p + m_e)` everywhere — assumes pure ionised hydrogen and ignores the inconsistency with the EOS ρ-definition
- **Where:** `src/opacity/magnetic_ff.jl:113, 154, 176` (`n_e = ρ / (m_p + m_e)`), `src/opacity/magnetic_modes.jl:30, 240` (`m_H = m_p + m_e`), and the EOS `ρ = m_p P/(2 k_B T)` everywhere.
- **What:** The EOS treats the gas as `n_e = n_p ≈ ρ/m_p` (since `2 (ρ/m_p) k_B T = P`). The opacity code treats `n_e = ρ/(m_p+m_e) ≈ ρ/m_p × (1 - m_e/m_p)`. The numerical difference is `m_e/m_p ≈ 5×10⁻⁴`, negligible. But it signals a missing single source of truth for "how do we convert mass density to electron number density?". When partial ionisation lands in Phase 4 (HANDOFF.md backlog), the two sites will need to be coupled and the discrepancy will matter.
- **Why it matters:** Latent bug class. The Project doc claims μ=0.5 — neither convention is exactly μ=0.5 in the m_e-correction sense.
- **Severity:** low (today), medium (Phase 4 partial ionisation).
- **Suggested fix:** Define `electron_density(ρ, T, x_H)` in a central place (PhysicalState module?) and call it everywhere.

### 23. `NSParams.f_col` exists but is dead in the new `render_spectral_cube` path
- **Where:** `src/pipeline/render.jl:35` (field); used at line 94 in the old `render_neutron_star` (modified-blackbody path) but ignored entirely in `render_spectral_cube`.
- **What:** A user setting `f_col=1.5` and rendering with the new pipeline gets the same answer as `f_col=1.0` — the colour correction is in the modified-blackbody placeholder, which the new pipeline replaces with grid lookup.
- **Why it matters:** Silently ignored parameters are footguns. Either remove from `NSParams` or warn.
- **Severity:** low.
- **Suggested fix:** Split `NSParams` into `StellarParams` (M, R, B_pole, T_pole, T_eq, obliquity, inclination, distance) and `BlackbodyRenderConfig { f_col, p }` used only by the legacy renderer.

### 24. `lookup_spectrum_scalar` is dead and slowly drifting from `lookup_spectrum`
- **Where:** `src/pipeline/atmosphere_grid.jl:144-163`. No callers.
- **What:** Copy of the per-frequency lookup. Will drift from the vector version as one is updated.
- **Severity:** low.
- **Suggested fix:** Delete.

### 25. Tests cannot exercise individual modules without an in-tree McPHAC checkout
- **Where:** `test/runtests.jl:74, 97`: `load_gaunt_table("refs/code/McPHAC/gffgu.dat")`.
- **What:** The Gaunt table loader has a hard dependency on a file inside a nested git repo (`refs/code/McPHAC`) which HANDOFF.md describes as "submodule-style, do not commit". Tests will fail in a CI environment that doesn't bring in McPHAC.
- **Why it matters:** No CI possible without setting up the McPHAC clone. Combined with the lack of artifacts (#14), the project cannot be built in a clean sandbox.
- **Severity:** medium (tests intersect with architecture here).
- **Suggested fix:** Vendor a 5KB subset of the Gaunt table into the package via Julia Artifacts, or include the file directly in `test/data/`.

---

## Things done WELL

1. **Module decomposition matches the physics pipeline.** EOS / TOV / surface / opacity / dielectric / atmosphere / RT / ray-tracing / colorimetry / pipeline is a clean, intuitive layering. A new reader can find any piece quickly.
2. **Physical constants are isolated, sourced (CODATA 2018), and re-exported once.** No constants are redefined in downstream modules. (Minor exception: `m_H = m_p + m_e` is defined in two files but the math is identical.)
3. **The B=0 limit recovers non-magnetic in the magnetic solver itself** (`magnetic_atmosphere.jl:289-297`). This is a real architectural choice (single code path with branch on B) rather than two duplicate solvers, even though most other dualities (`MagneticAtmosphereResult` vs `AtmosphereResult`, mode-aware vs mode-free `I_emergent`) lean the other way. Keep this discipline.
4. **Per-mode Planck splitting (B_ν/2) is consistently documented and applied.** The "Critical Physics Conventions" list in HANDOFF.md is reflected in code comments at every use site. This kind of cross-document discipline is rare and valuable.
5. **`vacuum_mode_parameters` separates concerns properly:** vacuum coefficients, polarization parameter, Kz derivation, and final mode weights are each their own small function. Easy to test piece-by-piece — which the test suite does.
6. **The adaptive Feautrier solver (`solve_feautrier_all_adaptive`) is implemented and packaged** even though it's not yet wired in. This is real "we found a better method" architecture-aware work, not a placeholder. Just needs to be made the default.

---

## Cross-cutting recommendation

Most of the highest-severity findings (#1, #2, #3, #7, #11) share a root cause: **the pipeline layer (`AtmosphereGrid`, `Renderer`) does not respect the abstractions of the physics layer.** The grid drops `θ_B` and mode information; the renderer doesn't know about gravity-from-(M,R); two functions ignore their own polarization output and average it away. Before Phase 4 lands a Kerr ray tracer that needs polarization along null geodesics, this layer needs a real `PhysicalSurfaceState → ObservedSpectrum` abstraction, parameterised over (T_eff, B, θ_B, [polarization]) and with explicit provenance and units. Designing that one interface well will pay off for the next decade of features.
