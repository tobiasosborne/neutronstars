# Code Review — NeutronStar.jl

Reviewer: Julia code smells and equation-vs-paper derivations.
Date: 2026-05-16, commit `95c76b3`.
Scope: Julia antipatterns, hot-loop allocations, magic numbers, derivations vs cited papers. Module architecture, test coverage, and tidiness are out of scope (handled by reviewers 01/02/04).

## Executive summary

- **Two real Coulomb-logarithm bugs** that the existing physics report (`docs/physics_report.tex`) missed: (a) `coulomb_magnetic.jl:102` falls back to `K0_val = 10.0` for small `x_n`, which is wrong by a factor ~2-3 (correct: `-log(x/2) - γ_E`); (b) the same file's commented "P&C 2003 Eq. (44a) integrates over `y ∈ [0, ∞)`" claim is now true at the integral level but `_Q_n_alpha` line 89 silently returns 0.0 when `log(base) * |n|` exceeds 300, which prematurely truncates the tail for marginally quantising fields.
- **`temp_correction.jl:197` recomputes `B̄_i` inside the inner k-loop** with a full O(K) sum, making the Rybicki temperature correction O(N·K²) instead of the O(N·K) it was designed to be — and the per-depth array `B_bar` precomputed on line 80-85 is *unused*. Equivalent magnetic code (`magnetic_atmosphere.jl:574`) does it correctly, so this is a clear regression in the non-magnetic path.
- **The Suleimanov flux correction `ΔT *= 1 + flux_damping × (flux_ratio^(-1/4) − 1)` is not from SPW09.** SPW09 §2 uses Avrett-Krook with their Eqs. (22-23) keyed on the per-depth energy-balance error `ε_Λ(m)`, plus a surface correction from Kurucz 1970. The code's grey-body scaling is heuristic; it is hardcoded as "Suleimanov-style" in the comments at `magnetic_atmosphere.jl:120` and in HANDOFF.md but no equation in SPW09 produces this formula. This is what produces the `[0.9, 1.1]` clamp on a quantity that should be a small correction.

The Phase 3b magnetic-RT formulae *that the report does cover* (P&C 2003 Eqs. 44b-e, 51, 52, 53 and the mode-vector setup) match the papers correctly. The real load-bearing issues are in the supporting infrastructure: numerical fallbacks, hot-loop allocations, and the gap between cited and implemented flux-correction schemes.

---

## Section A: Julia code smells

### A1. Heap-allocated 3-element `Vector{Float64}` in every opacity call
- **Where:** `src/opacity/dielectric_tensor.jl:277, 300, 320` (`return [wm1, w0, wp1]`, `return [0.0, 1.0, 0.0]`); consumed in `src/opacity/magnetic_modes.jl:242` (`w1, w2 = polarization_weights_vacuum(...)`).
- **What:** Both `compute_weights_from_K_Kz` and the special-case branches in `polarization_weights_*` return freshly heap-allocated `Vector{Float64}` of length 3. These are called from `_mode_cross_section_sum` in the hottest inner loop of the magnetic atmosphere — for a typical (N=200, K=50, 2 modes, 30 iterations) solve, that is ~600 000 calls, each allocating *two* 3-element vectors (one per mode).
- **Severity:** high (allocations dominate runtime in hot RT loops).
- **Suggested fix:** use `StaticArrays.SVector{3,Float64}` or return a `NTuple{3,Float64}`. The downstream consumer at `magnetic_modes.jl:246-248` already iterates with `enumerate((-1, 0, 1))` — replacing `w[idx]` with tuple indexing is one line.

### A2. `polarization_weights_vacuum` called twice per opacity sample
- **Where:** `src/opacity/magnetic_modes.jl:41-43, 57-59`; `src/atmosphere/magnetic_atmosphere.jl:281-282`.
- **What:** `mode_absorption(j, ν, θ_B, B, T, ρ)` and `mode_scattering(j, ν, θ_B, B, T, ρ)` are independent functions, but both call `_mode_cross_section_sum` which calls `polarization_weights_vacuum(ω, B, θ_B, n_e)` — an expensive function that internally also calls `stix_parameters` and `vacuum_mode_parameters`. In `_compute_magnetic_opacities!` these are called sequentially for each mode and depth-frequency point, so the polarization weights are evaluated **4× per (i, k, j)** instead of 1×.
- **Severity:** high (factor ~4 speedup available in the magnetic opacity path).
- **Suggested fix:** introduce `mode_opacity_split(j, ν, θ_B, B, T, ρ) → (κ_abs, κ_scat)` that returns both quantities from one weight call. Update `_compute_magnetic_opacities!` to use it.

### A3. `_nearest(grid, x) = argmin(abs.(grid .- x))` allocates per call
- **Where:** `src/pipeline/atmosphere_grid.jl:194-196`, called from `lookup_spectrum` line 127 and `lookup_spectrum_scalar` line 150.
- **What:** `abs.(grid .- x)` builds a brand-new `Vector{Float64}` of length `length(grid)` (≈50 for default config) on every lookup. A 256×256 render with 50 frequencies and bilinear (T,B) interpolation calls this 256·256·50·4 ≈ 13 million times.
- **Severity:** high.
- **Suggested fix:** the frequency grid is sorted; use `searchsortedfirst`/`searchsortedlast` and pick the closer of the two bracketing indices. Three lines, zero allocation.

### A4. `_feautrier_single` allocates `Vector{Matrix{Float64}}` per frequency, per mode
- **Where:** `src/atmosphere/magnetic_atmosphere.jl:381-383`; mirror in `src/atmosphere/feautrier.jl:218-220, 341-343`.
- **What:** `B_tilde = [zeros(M, M) for _ in 1:N_eff]`, plus `Q_tilde` and `C_store`. For (K=50, 2 modes) this is 100 calls per iteration, each allocating `3 × N_eff` matrices of size M×M, plus inner `copy(B)`/`copy(Q)`/`copy(C)` (lines 434, 438-439) and `B_tilde[i-1] \ Q_tilde[i-1]` (line 443) which allocates the result of `\\` twice (once unused, once used) per i. Profiling will show this dominates a magnetic solve.
- **Severity:** high.
- **Suggested fix:** preallocate `B_tilde`, `Q_tilde`, `C_store` once per `_solve_feautrier_mode!` call and reuse across frequencies. Use `lu!` (in-place LU) on a preallocated scratch matrix. The `copy(B)` calls (lines 434, 438) are unnecessary since `B` is freshly built every iteration.

### A5. `B̄_i` recomputed inside k-loop in non-magnetic Rybicki; precomputed `B_bar` array discarded
- **Where:** `src/atmosphere/temp_correction.jl:73-85` (precompute), `:197` (re-compute inline).
- **What:** Lines 73-85 precompute `B_bar[i] = Σ_k B_ν(ν_k, T_i) × κ_{i,k} × b_k / denom_i` once. Then line 197 *re-computes the same sum* inside `_build_rybicki_system!`, which is called `K` times per depth. Net cost: O(N·K²) instead of the intended O(N·K). The precomputed `B_bar` array is never read. The magnetic version (`magnetic_atmosphere.jl:574`) correctly uses the precomputed `B_bar[i]`.
- **Severity:** high (4-8× slowdown of the non-magnetic Rybicki step; K=50, N=100 means the per-iteration cost goes from 5000 sums to 250 000).
- **Suggested fix:** replace line 197 with `B_bar_i = B_bar[i]` and pass `B_bar` through the function signature like the magnetic version does.

### A6. `K0_val` fallback for tiny argument is numerically wrong
- **Where:** `src/opacity/coulomb_magnetic.jl:102`: `K0_val = x_n > 1e-10 ? besselk(0, x_n) : 10.0`.
- **What:** The comment says "K₀(x) ~ -ln(x) as x→0" — correct asymptotic — but the fallback returns 10.0 unconditionally. For `x_n = 1e-10` the true value is `-log(5e-11) - 0.5772 ≈ 23.1`, not 10. For `x_n = 1e-30` (the `θ → 0` longitudinal limit where the resonance argument vanishes) it should be ≈ 68. The 10.0 is more than 2× too small in the regime where the `α = ±1` contribution dominates.
- **Severity:** medium (only affects very low `(u - n β_e)`, which is near cyclotron resonance; but resonance points are exactly where one cares about accuracy).
- **Suggested fix:** replace the fallback with the correct asymptote `K0_val = -log(x_n / 2.0) - 0.5772156649` (Euler-Mascheroni). `besselk(0, x)` itself returns `Inf` for `x = 0`, so the guard is needed; but the fallback value must match the actual asymptote.

### A7. Inline `besselk(1, x)` fallback is also wrong
- **Where:** `src/opacity/coulomb_magnetic.jl:98`: `K1_val = x_n > 0 ? besselk(1, x_n) : 1.0`.
- **What:** The comment "K₁(x) ~ 1/x as x→0" is right, but the product `x_n × K₁(x_n)` (line 99) tends to 1.0 as `x_n → 0`, so returning `K1_val = 1.0` here makes `x_n × K1_val → 0` instead of → 1. The branch is reached only for `x_n == 0` exactly so it is rare, but the asymptotic limit is silently wrong by an O(1) factor for the entire `α = 0` longitudinal absorption channel at the cyclotron harmonic `u = n β_e`. Correct value: `x_n * K_1(x_n) → 1.0` so return `K1_val = 1.0 / x_n` if you want the product to be ≈ 1, or special-case the whole expression `A = 1.0 / (y + β_e/4)`.
- **Severity:** medium (same caveat as A6 — affects resonance points).
- **Suggested fix:** restructure `A_n^0 = x_n × K_1(x_n) / (y + β_e/4)` to special-case `x_n = 0` by setting `A = 1.0 / (y + β_e/4)` directly.

### A8. Premature truncation of Landau-level sum via `log(base) * |n|` check
- **Where:** `src/opacity/coulomb_magnetic.jl:86-90`: `if base <= 0 || log(base) * abs_n > 300.0; return 0.0; end`.
- **What:** For `β_e ≳ 1` and `y ≳ 1`, `sinh(β_e/2)` grows like `exp(β_e/2)/2`, so `base = (y+θ+ζ) × sinh(β_e/2)` can easily exceed `exp(5)` ≈ 150. For `|n| ≥ 60` then `log(150) × 60 > 300` and the integrand returns 0 even though the contribution should be `exp(-300)` ≈ 5e-131 — negligible numerically but the *sum* over `n` doesn't see that there's a smooth tail. More importantly, the truncation `N_max = min(50, ...)` on line 45 already caps |n| ≤ 50. So this guard catches nothing in the default path, but if `β_e` is small and `N_max = 50` the contributions from `|n| ≈ 50` *should* be exponentially small but the guard returns 0 abruptly rather than `exp(-large)`. Minor edge-case issue but the comment "Boltzmann suppression" is misleading — this is a *log-of-the-suppression* guard, not Boltzmann itself.
- **Severity:** low.
- **Suggested fix:** since the truncation is already by `N_max`, drop the `log(base) * |n| > 300` guard or replace it with `return exp(-log(base)*abs_n) * (A from the regular path)` so the integrand is smooth.

### A9. Inconsistent EOS formula between exported and private density update
- **Where:** `src/atmosphere/atm_structure.jl:144` (exported `update_atmosphere!`) vs `src/atmosphere/rt_atmosphere.jl:165` (private `_update_structure!`).
- **What:** Reviewer 01 already flagged this as an architecture-grade bug; from a code-smell angle the symptom is that the project has *two* implementations of "update density from (P, T) for ionised H" and they disagree by a factor of 2. Pattern is repeated again at `magnetic_atmosphere.jl:152` and `:200, :241`, each writing out `m_p * P[i] / (2.0 * k_B * T[i])` by hand. There is no `density_from_PT` helper in `AtmosphereStructure`. This is exactly the kind of duplicated formula that produces silent EOS bugs when one site is changed.
- **Severity:** high.
- **Suggested fix:** add `density_from_PT(P, T) = m_p * P / (2 * k_B * T)` to `AtmosphereStructure`, export it, replace all five inline copies.

### A10. `magnetic_modes.jl` and `magnetic_atmosphere.jl` each redefine `const m_H = m_p + m_e`
- **Where:** `src/opacity/magnetic_modes.jl:30`, `src/atmosphere/magnetic_atmosphere.jl:27`.
- **What:** Same constant, two modules. Not wrong, just a DRY violation. If the day comes when someone wants to switch to `m_u` or include neutron components, it'll be a hunt.
- **Severity:** low.
- **Suggested fix:** move to `PhysicalConstants` as `const m_H = m_p + m_e`.

### A11. Untyped return signatures on dielectric-tensor functions
- **Where:** `src/opacity/dielectric_tensor.jl:30, 66, 116, 129, 168, 236, 268, 297` — none declare `::Tuple{...}` or `::NTuple{...}` return types.
- **What:** `stix_parameters` returns `(S, D, P)`; `polarization_weights_*` returns `(Vector{Float64}, Vector{Float64})`; `vacuum_mode_parameters` returns a 6-tuple. With no return type annotation Julia infers them at compile time but inference can degrade if a code path conditionally returns differently-shaped objects (line 277/300/320 all return 3-element vectors; consistent, OK). The bigger issue is documentation: a caller has no way to know `vacuum_coefficients` returns `(a_hat, q, m)` without reading the source.
- **Severity:** low.
- **Suggested fix:** annotate return types `::NTuple{3,Float64}`, `::NTuple{6,Float64}`, etc. The compiler doesn't need it but humans and `@code_warntype` users do.

### A12. `lookup_spectrum` uses nearest-neighbour for frequency despite log-spaced grid
- **Where:** `src/pipeline/atmosphere_grid.jl:127, 150` (call to `_nearest`).
- **What:** Frequency interpolation is nearest-neighbour. The grid spans 0.05–120 `k_B T/h` over 50 points (log-spaced), so adjacent points differ by ~17% in frequency. The spectrum's peak (where the Wien displacement law puts most of the flux) varies by 17% between adjacent bins, and the nearest-neighbour quantization can produce visible step artifacts in `render_spectral_cube` outputs. The code does log-linear interpolation in `rt_emergent_spectrum` (line 228 in `rt_atmosphere.jl`); doing it in two places with two different methods is asking for trouble.
- **Severity:** medium (visible aliasing in spectral cubes, especially at low T_eff).
- **Suggested fix:** replicate the log-log interpolation pattern from `RTAtmosphere.rt_emergent_spectrum` in `lookup_spectrum`.

### A13. `flux_scale` clamp `[0.9, 1.1]` is a hack to mask divergence
- **Where:** `src/atmosphere/magnetic_atmosphere.jl:126`: `flux_scale = clamp(1.0 + flux_damping * (raw_flux_scale - 1.0), 0.9, 1.1)`.
- **What:** The clamp hides cases where the magnetic atmosphere is on a runaway path. With `flux_damping = 0.5` and ratio `0.1` (badly converging), `raw_flux_scale ≈ 1.78`, `1 + 0.5 × 0.78 = 1.39`, clamped to `1.1`. So the iteration silently caps the temperature correction at 10% even when the actual physics demands 40%. Combined with the temperature damper at line 132 (clamp to 30% max ΔT/T), this means the solver can declare "iteration limit reached" while still off the target by an order of magnitude in flux. The convergence test (line 139) requires `max_dT < tol AND |flux_ratio - 1| < flux_tol`, but if the damped flux scale always stays within 10%, the convergence flag depends solely on whether the *clamped* update happened to be small enough.
- **Severity:** medium-high (it works for `B = 10¹² G`, but the HANDOFF notes Suleimanov Fig 2 mismatch at `θ_B = 45°` for `B = 10¹⁴ G` — exactly the regime where this damper might be hiding a real iteration issue).
- **Suggested fix:** if `raw_flux_scale` is outside `[0.5, 2]`, log a warning and reduce `flux_damping` for the next iteration rather than silently clamping. Better: actually implement SPW09's Avrett-Krook correction (see B3).

### A14. `@assert` for hot-loop validation
- **Where:** `src/opacity/magnetic_modes.jl:236-237`, `src/opacity/magnetic_ff.jl:106-107`, `:171`, `src/opacity/dielectric_tensor.jl:31, 169`, `src/opacity/hydrogen_ff.jl:28`, etc.
- **What:** Asserts run on every call in tight inner loops. With Julia's default `--check-bounds=yes` they cannot be elided. The `@assert ν > 0 && B >= 0 && T > 0 && ρ > 0` style triggers comparisons + short-circuit evaluation per call — measurable in inner loops.
- **Severity:** low.
- **Suggested fix:** convert input-domain checks to one-time validation in the *outer* solver loop. Keep `@assert` only for invariants you genuinely want stripped at `-O3`.

### A15. `@printf` peppered through iterations; no logging level
- **Where:** `src/atmosphere/rt_atmosphere.jl` (10 sites), `src/atmosphere/magnetic_atmosphere.jl` (6 sites), `src/pipeline/atmosphere_grid.jl` (5 sites), `src/pipeline/render.jl` (8 sites).
- **What:** Everything is `verbose && @printf(...)`. A solve produces ~60 lines of stdout that look fine in REPL but bury anything important. No `@info`/`@warn`/`@debug` levels. The single `@warn` in the codebase is in `tov.jl:221`. No way to silence a single iteration's per-frequency diagnostics while keeping the per-iteration summary.
- **Severity:** medium (user pain, not physics).
- **Suggested fix:** replace `verbose &&` with `@info` / `@debug` and let users set the log level via `Logging.global_logger`. Or thread a `LogLevel` enum through the call signatures.

### A16. Doc-comment for `polarization_weights_full` is identical to `polarization_weights_cold`'s and labels the wrong function
- **Where:** `src/opacity/dielectric_tensor.jl:50-65` (the docstring above `polarization_weights_cold`) and `:112-114` (the docstring above `polarization_weights_full`).
- **What:** The docstring at line 50 says "`polarization_weights_full(ω, B, θ_B, n_e)`" but is attached to `polarization_weights_cold`. The `polarization_weights_full` function (line 116) is now just a compatibility shim. This is a stale rename: the docstrings were not updated when the function was split.
- **Severity:** low.
- **Suggested fix:** rewrite the docstring at line 50 to describe `polarization_weights_cold` and reference the q-based derivation; update the line 112 docstring to clearly mark `polarization_weights_full` as DEPRECATED.

### A17. Cross-module name collision: `K` is both an integer (frequency count) and a polarization ratio
- **Where:** Throughout `magnetic_atmosphere.jl` (`K = length(ν_grid)`), `dielectric_tensor.jl` (`K1, K2, Kz`), `feautrier.jl` (`K_k`, `K_k_vector`).
- **What:** A reader/agent encountering `K` in `_build_rybicki_system_magnetic!` has to mentally disambiguate: is this the loop-count `K`, the polarization ratio `K_j`, or the Rybicki RHS vector `K_k`? Same letter, three meanings.
- **Severity:** low.
- **Suggested fix:** rename `K` (frequency count) to `nν` everywhere; rename the Rybicki RHS to `rhs_k` or `Q_k`.

### A18. `argmin(abs.(grid .- x))` pattern is the wrong stdlib reach
- **Where:** `src/pipeline/atmosphere_grid.jl:195`.
- **What:** Same issue as A3 but worth mentioning as a separate Julia antipattern. Whenever you see `argmin(abs.(...))` in performance-critical code, the right tool is `searchsorted*` on the (sorted) grid. This is a textbook "reinventing stdlib" smell.
- **Severity:** medium (covered by A3 fix).
- **Suggested fix:** see A3.

### A19. `Float64` hardcoded throughout — no AD or `Float32` path
- **Where:** Every function signature in `opacity/`, `atmosphere/`, `eos/`. Examples: `kappa_ff(ν::Float64, T::Float64, ρ::Float64, ...)::Float64`, `polarization_weights_vacuum(ω::Float64, ...)`.
- **What:** All numeric arguments and return types are concretely `Float64`. If a user wanted to differentiate through the atmosphere solver (ForwardDiff, Enzyme, Zygote) to e.g. compute `dT_eff/d(observed_spectrum)` or to fit observed neutron-star data, they'd have to refactor every signature. Likewise no `Float32` path for memory-constrained renders.
- **Severity:** low-medium (depends on whether AD is on the roadmap).
- **Suggested fix:** parameterize on `T<:Real`: `kappa_ff(ν::T, T_temp::T, ρ::T, gaunt)::T where {T<:Real}`. Most of the code already does `2π * ν`, `ω^2`, etc. which are type-stable polymorphic. The hot-loop functions in `magnetic_ff.jl` and `dielectric_tensor.jl` are the main beneficiaries.

### A20. `eachindex(μ_full)` followed by index-based access of `μ_full[i]`, `w_full[i]`
- **Where:** `src/atmosphere/feautrier.jl:32-37` (`gauss_legendre_half`).
- **What:** Iterating with `eachindex(μ_full)` then pushing into `μ` and `w` is fine semantically but the function's hot path is at solver setup, not per-iteration, so this is forgivable. More concerning is that `_gauss_legendre` (line 457-463) builds a `SymTridiagonal`, calls `eigen`, and reads `vecs[1, :].^2` — this is the Golub-Welsch algorithm but does an `eigen` call (which uses LAPACK and returns sorted eigenvalues). For M=8, this is overkill — there are precomputed Gauss-Legendre weights for small M in `FastGaussQuadrature.jl`. But more importantly, this is called once per `solve_atmosphere` — not a perf issue.
- **Severity:** low.
- **Suggested fix:** add a small lookup table for M ∈ {4, 6, 8, 10, 12, 16} (the practical range), or depend on `FastGaussQuadrature.jl`.

### A21. `mutable struct AtmosphereColumn` with all abstract field types
- **Where:** `src/atmosphere/atm_structure.jl:26-43`.
- **What:** The struct is mutable and uses concrete field types (`Vector{Float64}`, `Matrix{Float64}`) — that part is fine. But `mutable struct` adds heap allocation per instance even though the column itself is built once per solve. The mutation pattern (`col.T .= T_new`) only mutates the *contents* of the contained arrays, never reassigns the array references. So `mutable` is unnecessary.
- **Severity:** low.
- **Suggested fix:** make it `struct` (immutable). The `.= T_new` pattern still works on `Vector{Float64}` fields. The only place that actually reassigns a field is `update_atmosphere!` line 142 (`col.T .= T_new`) which is `.=` (broadcast into existing array), so immutable is fine.

### A22. `m_H = m_p + m_e` used as proton mass in `n_e = ρ / m_H`
- **Where:** `src/opacity/magnetic_modes.jl:240`; `src/opacity/magnetic_ff.jl:113, 154, 176, 183`.
- **What:** For fully ionised H, `n_e = n_p = ρ/(m_p + m_e)` is correct (the *total* mass per electron-proton pair is `m_H = m_p + m_e`). This is the convention used in H12 and PC03. The code is right but the variable name `n_e = ρ / m_H` reads like a per-atom density and a reader unfamiliar with the EOS convention might think `m_H` should be `m_p` (which would give a 0.05% high density). A one-line comment would prevent confusion.
- **Severity:** low.
- **Suggested fix:** add `# n_e = n_p for fully ionized H; m_H = m_p + m_e is the per-pair mass` at first use.

### A23. `besselk(0, u/2)` overflow guard is too tight
- **Where:** `src/opacity/coulomb_magnetic.jl:120-124` and `src/opacity/magnetic_ff.jl:38-42`.
- **What:** Both files guard `K_0(x)` evaluation with `if u/2 > 300; return asymptote`. For `u = h ν / (k_B T)` at `T = 5×10⁵ K` and `ν = 10²² Hz` (hard X-ray well beyond NS spectrum), `u ≈ 1900`, so `u/2 ≈ 950 > 300`. The asymptote `√(π/2x)` is correct, but `magnetic_ff.jl:39` uses `if u2 > 500; return 1.0` which is the *wrong asymptote* (should be `√(π/(2 × 250))` ≈ 0.079, not 1.0). Function `coulomb_log_classical` (the non-safe version) at line 37-43 of `magnetic_ff.jl` therefore returns 1.0 at very high u, while `coulomb_log_classical_safe` in `coulomb_magnetic.jl:118-125` returns the correct asymptote. These are two implementations of the same equation that disagree at the asymptotic limit.
- **Severity:** medium (affects very high frequency or low T, but is silently wrong).
- **Suggested fix:** delete `coulomb_log_classical` (line 37 of `magnetic_ff.jl`) and have all callers use `coulomb_log_classical_safe`. There's no reason to keep two versions.

### A24. `_log_interp` clamps to log10 of 1e-30 — silent zero clipping
- **Where:** `src/atmosphere/feautrier.jl:302-313`.
- **What:** `logx_old = log10.(max.(x_old, 1e-30))` — if any opacity or albedo is less than 1e-30 (possible for `ρ_alb` at the surface where scattering can dominate by huge factors), it's clamped silently. Combined with the linear interpolation, this means a deep-grid value of, say, `ρ_alb = 1e-50` gets treated as `1e-30`, which can pollute the per-frequency Feautrier solve. No warning. The downstream `clamp(f_common[i], 0.0, 1.0)` on line 191 hides the symptom.
- **Severity:** low.
- **Suggested fix:** make the clamp value an argument with a sensible default (say `eps(Float64)`) and assert the input is non-negative.

### A25. `_solve_single_frequency` has dead code branches and confused bottom-BC logic
- **Where:** `src/atmosphere/feautrier.jl:363-391`.
- **What:** Lines 370-373 set `B_tilde[i] = B; Q_tilde[i] = Q`. Then lines 372-380 conditionally compute `A_mat` but never use it ("Actually at bottom BC: P = B_ν regardless, so just set it"). Then lines 382-383 set `B_tilde[i] = B; Q_tilde[i] = Q` *again*. Lines 386-390 then run *another* `if i > 1` branch that sets `B_tilde[i] = B; Q_tilde[i] = Q` for the *third* time, with the comment "Keep as identity / Keep as B_ν". The code is salvageable — the LTE bottom BC `P_N = B_ν` is correct — but the dead branches and recursive overwrites are a sign of incomplete refactoring. Compare with the cleaner magnetic version at `magnetic_atmosphere.jl:402-411`.
- **Severity:** low (code is correct, just confusing).
- **Suggested fix:** delete lines 372-391 entirely; just keep `for j in 1:M; B[j,j] = 1.0; end; Q .= Bν; B_tilde[i] = B; Q_tilde[i] = Q; C_store[i] = C; continue`.

---

## Section B: Derivation discrepancies

### B1. `_magnetic_eddington_temperature` uses non-magnetic surface BC `T(0) = 0.265 T_eff`
- **Where:** `src/atmosphere/magnetic_atmosphere.jl:224`: `T[1] = 0.265 * T_eff`.
- **Code's formula:** Surface temperature initialised to `0.265 T_eff`, matching the McPHAC convention for `B = 0`.
- **Paper formula:** SPW09 §2 (paragraph 4 of "Scheme of calculations") starts the iteration from "a grey temperature distribution" but does not specify a magnetic surface BC. For `B → ∞`, the surface limit of the grey atmosphere is *not* `0.265 T_eff`: the grey factor is `(3/4)^{1/4}(τ + 2/3)^{1/4}`, evaluated at `τ_R = 0`, gives `T = (1/2)^{1/4} T_eff = 0.84 T_eff` for a pure-absorption atmosphere. The McPHAC value `0.265` (from H12 §3.1) is empirical — derived by integrating their non-magnetic Eddington equation `dT/dy = (3/16) k_R (T/T_eff)³` down to the surface and reading off the result. The same integration with the *magnetic* Rosseland `K_⊥` would give a different surface limit because `K_⊥` ≠ `k_R` near the surface where ionisation is partial. Using `0.265` for both is an unjustified import of the non-magnetic answer.
- **Match:** the surface BC is wrong by ~3× in the worst case for a strongly magnetised atmosphere.
- **Severity:** medium (affects only the initial guess; the Rybicki iteration may recover the right surface T, but the initial guess influences convergence speed and could explain the `B = 10¹⁴`, `θ_B = 45°` slow convergence noted in HANDOFF.md).
- **Fix:** integrate the magnetic Eddington equation `T⁴ = (3/4) T_eff⁴ (τ + 2/3)` from `τ = ∞` down to `τ = 0` *self-consistently* with the magnetic Rosseland mean, rather than hard-coding `0.265`. Alternatively, set `T[1] = (1/2)^{1/4} T_eff = 0.841 T_eff` as the analytic grey-atmosphere surface limit and let the iteration drift to the right value.

### B2. Mode-vector formula in `compute_weights_from_K_Kz` — sign convention unverified
- **Where:** `src/opacity/dielectric_tensor.jl:297-321`.
- **Code's formula:** Given (K, Kz), with `K = E_⊥2 / E_⊥1` and `Kz = E_z / E_⊥1`,
  ```
  transverse = K cosθ + Kz sinθ
  longitudinal = K sinθ - Kz cosθ
  |e_+1|² = (transverse + 1)² / (2N)
  |e_-1|² = (transverse - 1)² / (2N)
  |e_0|²  = longitudinal² / N
  N = 1 + K² + Kz²
  ```
- **Paper formula:** Per HANDOFF.md the authoritative source is "van Adelsberg & Lai 2006 B-frame cyclic component formulae" but the local `refs/ho_lai_2001_vacuum_pol.pdf` is mislabelled. P&C 2003 Eq. 26 gives the cyclic weights via `(K_j ± 1)² / (2(1 + K_j² + K_{z,j}²))` for the `± 1` components and `K_{z,j}² / (1 + K_j² + K_{z,j}²)` for the 0 component. The code's `transverse` and `longitudinal` variables suggest a rotation from the wave-vector frame to the B-frame, which is the van Adelsberg & Lai (2006) Eq. (16)-(18) construction. *I cannot verify this against the paper because the referenced PDF is wrong.*
- **Match:** unverified.
- **Severity:** **NEEDS PAPER VERIFICATION**. This is the formula at the heart of the unresolved `θ_B = 45°` mismatch (HANDOFF.md "Immediate Next Task"). The HANDOFF.md regression test (`K1 = -0.025737`, `K2 = 38.849`, `K1 × K2 ≈ -1`) checks the *intermediate* `K` values but not the final `|e_α|²` against a paper-derived reference. Suggest the next step explicitly: download van Adelsberg & Lai (2006) from arXiv, re-derive `compute_weights_from_K_Kz` from their Eqs. (16)-(18), and add a regression test for the cyclic weights at three benchmark `(B, ω, θ_B, ρ)` points pulled from their Fig. 1 or Fig. 2.

### B3. The flux correction is not Suleimanov-style
- **Where:** `src/atmosphere/magnetic_atmosphere.jl:119-128`.
- **Code's formula:**
  ```
  F_bol = _bolometric_flux_2mode(P_all, μ, w, ν_grid)
  flux_ratio = F_bol / (σ_SB × T_eff^4)
  raw_flux_scale = flux_ratio^(-0.25)
  flux_scale = clamp(1 + flux_damping × (raw_flux_scale - 1), 0.9, 1.1)
  ΔT .+= (flux_scale - 1) × T
  ```
- **Paper formula:** SPW09 §2, third paragraph after Eq. (21): "Temperature corrections are then evaluated using three different procedures. The first is the integral Λ-iteration method, modified for the two-mode radiation transfer, based on the energy balance equation (9). In this method the temperature correction for a particular depth is found from `ΔT_Λ(m) = -ε_Λ(m) / ∫Φ_ν dν`. The second procedure is the Avrett-Krook flux correction, which uses the relative flux error `ε_H(m)` ... And the third one is the surface correction, which is based on the emergent flux error". SPW09 Eq. (22) gives the Λ-iteration correction. SPW09 Eq. (20) defines `ε_H(m) = 1 - H_0 / ∫₀^∞ H_ν(m) dν`. The Avrett-Krook formula is *spatially resolved* per-depth using the deep flux error.
- **Match:** The code's formula is *globally* averaged (uses only the surface `F_bol`) and applied uniformly to all depths via `(flux_scale - 1) × T`. SPW09's Avrett-Krook uses the depth-dependent `ε_H(m)` to construct a depth-dependent correction. The clamp to `[0.9, 1.1]` and the `flux_damping = 0.5` factor are heuristic stabilisers with no SPW09 analogue. The exponent `-1/4` is a heuristic grey-body scaling (motivated by `T_eff ∝ F^(1/4)`) but is not the prescription anywhere in SPW09.
- **Severity:** high — the code claims to implement SPW09 in HANDOFF.md ("Damped Suleimanov-style flux correction") and in the inline comment at line 119-121, but it does not.
- **Fix:** either (a) rename "Suleimanov-style" to "ad-hoc grey scaling" everywhere and document the heuristic origin, or (b) implement the actual Avrett-Krook scheme: compute `H_ν(m_i)` at every depth, form `ε_H(m_i)`, and apply the temperature correction `ΔT_AK(m_i) = α × ε_H(m_i) × T_i / (∂ ln B_ν/∂ ln T)` per Kurucz 1970 / SPW09. The Suleimanov Fig 2 `θ_B = 45°` mismatch may partly be explained by this — the ad-hoc clamp prevents the temperature profile from adjusting enough at depth to balance the mode-asymmetric opacity gradient.

### B4. Magnetic Feautrier scattering kernel is isotropic, not anisotropic dipole
- **Where:** `src/atmosphere/magnetic_atmosphere.jl:428-430` vs `src/atmosphere/feautrier.jl:411-419`.
- **Code's formula (magnetic):** `B[j, jp] -= 2.0 * ρ * 0.5 * w[jp]` — isotropic kernel `dσ/dμ' = 1/2`.
- **Code's formula (non-magnetic):** `B[j, jp] -= 2.0 * ρ_alb * (3/16)(3 + 3μ²μ'² - μ² - μ'²) * w[jp]` — Hummer (1962) dipole phase function.
- **Paper formula:** SPW09 Eq. (8) describes `σ_ν^{i,j}(μ, φ; μ', φ')` as a full *inter-mode and anisotropic* scattering matrix. The natural fully magnetic kernel for Thomson scattering on free electrons is the magnetised Compton kernel (e.g. Mészáros 1992), which for low-frequency limit reduces to the dipole pattern weighted by polarisation projections.
- **Match:** Code uses *isotropic* for magnetic; the non-magnetic code uses *anisotropic dipole*. This is acknowledged in HANDOFF.md "Still Not Physically Complete → Scalar / local scattering approximation". But the inconsistency means: in the B → 0 limit, the magnetic solver does *not* reduce to the non-magnetic solver because the scattering kernels differ. The HANDOFF claim "B=0 limit recovers non-magnetic to <0.1%" is suspicious — if true, it's only because the scattering contribution is small at the test conditions, not because the kernels agree.
- **Severity:** medium (documented as a known limitation; severity comes from the inconsistency in the B→0 limit which contradicts a verification claim).
- **Fix:** at minimum, make the magnetic Feautrier use the same Hummer dipole kernel as the non-magnetic, weighted by `ρ_alb`. Long-term: implement the SPW09 Eq. (8) inter-mode redistribution.

### B5. `_bolometric_flux_2mode` uses trapezoid `dν = ν[k+1] - ν[k]` but the per-mode `I_emergent = 2P_surface` factor pulls in a `2 × 2π × ...` weighting that may double-count
- **Where:** `src/atmosphere/magnetic_atmosphere.jl:591-603`.
- **Code's formula:**
  ```
  F = Σ_k Σ_j Σ_mode 2π × μ_j × 2.0 × P_all[1,j,k,mode] × w_j × dν
  ```
- **Paper formula:** SPW09 Eq. (18): `4π ∫ H_ν(0) dν = σ T_eff⁴`, with Eq. (19): `H_ν(m) = (1/4π) Σ_{i=1,2} ∫ dφ ∫ μ I_ν^i dμ`. So:
  ```
  σ T_eff⁴ = 4π × (1/4π) Σ_{i=1,2} ∫ ∫ ∫ μ I_ν^i dμ dφ dν
           = Σ_{i=1,2} 2π × ∫ dν ∫ μ I_ν^i dμ          (after φ-integration over 2π)
           = Σ_{i=1,2} 2π × ∫ dν × (sum_j μ_j a_j I_emergent_{j,k,i})
  ```
  With `I_emergent = 2 P_all[1,j,k,mode]` (from the Feautrier `I = 2P` convention), this becomes
  ```
  σ T_eff⁴ = Σ_{i=1,2} 2π × Σ_k dν × Σ_j μ_j w_j × 2 × P_all[1,j,k,mode]
  ```
  The code matches.
- **Match:** ✓ correct.
- **Severity:** none. Listed here only because the same expression in HANDOFF.md ("F = 2π ∫₀¹ I(μ) μ dμ (not 4π)") was potentially ambiguous about whether to sum over modes outside or inside the 2π factor. The code is right.

### B6. `cold_polarization_parameter` falls back to `copysign(Inf, ...)` instead of using the limiting analytic form
- **Where:** `src/opacity/dielectric_tensor.jl:196-200`.
- **Code's formula:**
  ```
  if abs(D) < 1e-30 || abs(cosθ) < 1e-30
      return copysign(Inf, (P - S) * D * cosθ)
  end
  return (P - S) * sinθ^2 / (2 * D * cosθ)
  ```
- **Paper formula:** P&C 2003 Eq. 25, exact: `q = (P - S) sin²θ / (2 D cosθ)`. For `θ → π/2` (cosθ → 0), `q → ∞`, and the modes become `K₁ → 0` (transverse → 0) and `K₂ → -∞`, which in the cyclic weights gives `|e_{1,0}|² → 1` (longitudinal) and `|e_{2,±1}|² → 1/2` each (transverse circular mix). The code handles this special case at `polarization_weights_cold:81-85` before calling `cold_polarization_parameter`, so the `Inf` return is never reached *in production*. But `vacuum_polarization_parameter` (line 220) calls `cold_polarization_parameter` as a fallback — and if a caller passes `cosθ = 1e-31` (just below the special-case threshold of 1e-8 in `polarization_weights_cold` but in the regime where `1e-30` triggers), they get `Inf` propagated through `vacuum_mode_parameters`. That `Inf` will cascade through `vacuum_mode_parameters → vacuum_mode_parameters → compute_weights_from_K_Kz`, which has its own guard `!isfinite(K) || !isfinite(Kz)` returning `[0.0, 1.0, 0.0]` on line 299-301. So the final answer is correct *by accident*. Smelly defensive programming.
- **Match:** correct numerically.
- **Severity:** low (correct by accident).
- **Fix:** make `cold_polarization_parameter` either error or return a `Union{Float64, Missing}` for ill-defined inputs; let callers decide on the limiting behaviour.

### B7. `vacuum_coefficients` exponent on `b^(3/4)` and `b^(5/4)` not double-checked against Potekhin, Lai & Chabrier (2004)
- **Where:** `src/opacity/dielectric_tensor.jl:175-185`.
- **Code's formula:**
  ```
  a_hat = -(2 α / 9π) × log(1 + (b²/5) × (1 + 0.25487 b^(3/4)) / (1 + 0.75 b^(5/4)))
  q     = (7 α / 45π) × b² × (1 + 1.2 b) / (1 + 1.33 b + 0.56 b²)
  m     = -(α / 3π) × b² / (3.75 + 2.7 b^(5/4) + b²)
  ```
- **Paper formula:** Potekhin, Lai & Chabrier (2004) "Electromagnetic Polarisation in Partially Ionised Plasmas..." Appendix A, Eqs. (A7)-(A9). The local PDF `refs/potekhin_chabrier_ho_2014_opacities.pdf` is *mislabelled* (it's actually Boev & Kovalev 2014 per HANDOFF.md). The HANDOFF.md regression test at `B = 10¹⁴ G` gives `a_hat = -2.0563758...×10⁻⁴`, `q_vac = 1.0012983×10⁻³`, `m_vac = -2.4251×10⁻⁴`. These values are *self-consistent with the code* but I cannot cross-reference against the actual paper.
- **Match:** **NEEDS PAPER VERIFICATION**. Reproduce the figures of PLC2004 from the arXiv source `/tmp/potekhin_lai_chabrier_2004_src/kke.tex` (per HANDOFF.md) and add a regression test against their Fig. 1 or a tabulated value.
- **Severity:** medium (load-bearing for all vacuum-polarisation effects; if the constants `0.25487`, `0.75`, `1.2`, `1.33`, `0.56`, `3.75`, `2.7` are off in even one digit, the X-mode opacity near the vacuum resonance is wrong).

### B8. Magnetic Coulomb log integral upper bound is `Inf` but uses `maxevals = 2000`
- **Where:** `src/opacity/coulomb_magnetic.jl:54`: `val, _ = quadgk(integrand, 0.0, Inf; rtol=1e-4, maxevals=2000)`.
- **Code's formula:** P&C 2003 Eq. 44a, `Λ_α^ff = (3/4) e^(u/2) Σ_n ∫_0^∞ Q_n^α dy`. The code passes `Inf` as the upper limit (good — earlier `y = 20` cutoff was the documented Eq. 44e bug).
- **Paper formula:** P&C 2003 Eq. 44a explicitly: integral over `y ∈ [0, ∞)`. The longitudinal channel (α = 0) for low-energy photons has a long power-law-like tail in `y` that the previous `y = 20` cutoff missed.
- **Match:** ✓ correct in intent. But: `quadgk` with `maxevals = 2000` and a semi-infinite domain may not converge to the requested `rtol = 1e-4` for sharply-peaked integrands near `x_n` resonances. The code does not check `_, err = quadgk(...)`; the discarded error estimate could be large. No `verbose` output if `quadgk` hits `maxevals`.
- **Severity:** medium.
- **Fix:** capture the error estimate: `val, err = quadgk(...)`; if `err / max(abs(val), 1e-10) > 1e-3`, log a warning. Also, for `α = 0` consider splitting the integral as `∫_0^1 + ∫_1^∞` to give `quadgk` better adaptive behavior.

### B9. `mode_absorption + mode_scattering ≠ mode_opacity` may fail in vacuum-resonance regime
- **Where:** `src/opacity/magnetic_modes.jl:39-72`.
- **Code's formula:** `mode_opacity = mode_absorption + mode_scattering`, both summing `|e_{j,α}|²` weights from `polarization_weights_vacuum`.
- **Paper formula:** P&C 2003 Eq. 27 for the *total* mode opacity uses `σ_α^tot = σ_α^ff + σ_α^pp + σ_α^scat`. The decomposition into absorption + scattering is the code's own contribution; P&C 2003 don't write it separately. The split is correct *if* `|e_{j,α}|²` is the same across all three cross-section types. It is — there's no polarization-dependent cross-section weighting.
- **Match:** ✓ correct.
- **Severity:** none. The test in `test/runtests.jl:20` `@test κ_total ≈ κ_abs + κ_scat rtol=1e-12` verifies this identity.

### B10. Mass/density used to compute `n_e` differs between `nu_ff_alpha` and `_sigma_ff_nonmag`
- **Where:** `src/opacity/magnetic_ff.jl:113` (`n_e = ρ / (m_p + m_e)`) vs `:154` (`n_e = ρ / (m_p + m_e)`) — actually both use the same thing.
- Looking again: both use `ρ / (m_p + m_e)`. ✓ No issue. (I was checking for the `m_p` vs `m_H` mistake noted in A22.)
- **Match:** ✓ correct.
- **Severity:** none.

### B11. `σ_T / 3` non-magnetic fallback in `sigma_scat_alpha` is correct only after weight integration
- **Where:** `src/opacity/magnetic_ff.jl:197-200`.
- **Code's formula:** `if B < 1e6; return σ_T / 3.0; end` for each polarisation α.
- **Paper formula:** P&C 2003 Eq. 33 in the B→0 limit gives `σ_α^{s,e} → σ_T` (all polarisations) since `(ω + α ω_ce)² → ω²`. So per-polarisation cross-section is `σ_T`, *not* `σ_T/3`.
- **Match:** **discrepancy** at first glance. *But* the code's `_mode_cross_section_sum` does `κ_mode = Σ_α |e_α|² × σ_α`. Sum of `|e_α|²` over α = 1. So `κ_mode = σ_α^scat` per polarisation (since all are equal at B=0) = `σ_T`. With the `σ_T/3` fallback, this becomes `σ_T/3` per mode → mode opacity is `σ_T/3` per mode → effective opacity (harmonic mean) is `σ_T/3`. That's *3× too low*.
- **Match:** the `σ_T/3` is a holdover from an earlier convention where the three polarisations were summed (giving `3 × σ_T/3 = σ_T` total). Under the current `Σ_α |e_α|² × σ_α` scheme with `Σ_α |e_α|² = 1`, the fallback should be `σ_T` per polarisation, not `σ_T/3`.
- **Severity:** medium-high. This branch is only triggered for `B < 1e6 G`, which is below any neutron-star scenario. So in practice it's dead code. But it *is* dead code that gives wrong answers, and `mode_scattering` is called with `B = 0` for the non-magnetic-comparison branch in `_compute_magnetic_opacities!:289-298` — except that branch *doesn't* call `sigma_scat_alpha` (it uses `sigma_thomson()` from the non-magnetic module). So in practice the bad branch is not exercised. Still worth fixing.
- **Fix:** change `return σ_T / 3.0` to `return σ_T`, document why (the weight sum is now `1.0`, not `3.0`). Add a regression test: for `B → 1e6` (just above threshold), the effective opacity should match `σ_T / (m_p + m_e)` to within Doppler broadening.

### B12. Stix `D` sign convention differs from textbook
- **Where:** `src/opacity/dielectric_tensor.jl:43`.
- **Code's formula:** `D = ω_ce × ω_pe² / (ω × de) - ω_cp × ω_pi² / (ω × dp)` where `de = ω² - ω_ce²`, `dp = ω² - ω_cp²`.
- **Paper formula:** Stix (1992) "Waves in Plasmas" Ch. 1, standard form: `D = Σ_s (q_s/|q_s|) × ω_ps² × ω_cs / [ω(ω² - ω_cs²)]`. For electrons (`q = -e`, sign factor `-1`) the contribution is `-ω_ce ω_pe²/(ω(ω²-ω_ce²))`. For protons (`q = +e`, sign factor `+1`) it is `+ω_cp ω_pi²/(ω(ω²-ω_cp²))`. Sum: `D_textbook = -ω_ce ω_pe²/(ω×de) + ω_cp ω_pi²/(ω×dp) = -D_code`.
- **Match:** Code's `D` is the *negative* of the textbook Stix `D`. The downstream `q = (P-S) sin²θ / (2 D cosθ)` then has opposite sign of textbook. The mode quadratic `K² + 2qK - 1 = 0` has roots `K_{1,2} = -q ± √(1+q²)`, so flipping the sign of `q` swaps `K_1 ↔ -K_2` (and the labels of "extraordinary" and "ordinary"). The cyclic weights `(K ± 1)²/(2N)` are *not* invariant under `K → -K` (they swap `α = +1 ↔ -1`).
- **Severity:** medium — but: P&C 2003 follow Ginzburg (1970) conventions, which differ from Stix. Without comparing P&C's actual sign of `D` to the code I cannot say if this is a bug or just a different sign convention. The fact that the cold-plasma tests pass and that the B=0 limit recovers Thomson suggests the convention is internally consistent. **NEEDS PAPER VERIFICATION**: locate Ginzburg (1970) §10 or P&C 2003 Eq. 22-25 in their source and confirm the sign of `D`. Cross-reference the code's `q` against their `q̄(ω) sin²θ/(2 cosθ)` expression of Eq. 26 — the sign should match.
- **Fix:** if the convention is intentional, add a one-line comment `# Sign of D follows P&C 2003 / Ginzburg, opposite of Stix textbook`. If it's a bug, flip the sign and watch all the regression tests.

### B13. Greenstein-Hartke `T(θ) = T_eq + (T_pole - T_eq) cos²θ` — paper says different
- **Where:** `src/surface/dipole.jl:42`.
- **Code's formula:** `T = T_eq + (T_pole - T_eq) cos²(θ_B)`.
- **Paper formula:** Greenstein & Hartke (1983) "Pulselike character of blackbody radiation from neutron stars" derives the temperature distribution for two-pole NS heating. Their Eq. (10) / (11) gives the surface temperature as `T(θ) = T_p × |cos θ|^{1/4}` (where `T_p` is the polar temperature) for a pure dipole field with anisotropic thermal conductivity. The code's `T_eq + (T_pole - T_eq) cos²θ` is an *interpolation form*, not the GH 1983 result.
- **Match:** **discrepancy**. The code is a phenomenological cos²(θ) interpolation between two parameters; GH 1983's actual prediction is `T(θ) ∝ |cos θ|^{1/4}` (a *much* shallower distribution).
- **Severity:** medium (this is the surface temperature map used in every render; the choice between `cos²` and `cos^{1/4}` changes the equatorial-to-polar temperature ratio).
- **Fix:** either rename the function and remove the GH attribution (it's just a `T_eq + ΔT cos²θ` ansatz), or implement the actual GH 1983 result. The current docstring says "Source: Greenstein & Hartke (1983) Eq. 1" which is not correct — their Eq. 1 is the dipole field magnitude, not the temperature.

---

## Things done well

1. **Magnetic Coulomb logarithm integration domain fix.** Switching from `[0, 20]` to `[0, ∞)` for the `quadgk` upper limit in `coulomb_magnetic.jl:54` is the right fix for the longitudinal polarisation tail. Documented in HANDOFF.md and the inline comment.
2. **Mode-separated angle averaging.** `magnetic_modes.jl:102-149` correctly angle-averages each mode independently before harmonic-mean over modes, matching P&C 2003 Eq. 30. The header comment explicitly calls out the previous bug ("Previous code combined modes first (wrong)").
3. **`vacuum_coefficients` returns analytic fits to PLC2004 A7-A9** with regression tests at the canonical `B = 10¹⁴ G` test point (`test/runtests.jl:46-48`). The intermediate `(β, r, K1, K2, Kz1, Kz2)` are also pinned. This is *exactly* the right way to lock in a paper-derived numerical recipe.
4. **Surface BC discretization in `temp_correction.jl:207`** uses the McPHAC `½×(k₁+k₂)` convention even though the interior `Δ` formula uses `(k_i + k_{i-1})` without the half. The comment explicitly notes the inconsistency between H12 Eq. A12 and A28-A30. Compatible with McPHAC source.
5. **Eq. 44e fix and inline regression test.** `coulomb_magnetic.jl:74` implements `x_n = |u - n β_e| × √(1/4 + y/β_e)` correctly; `test/runtests.jl:24-36` pins `Λ_0 ≈ 7.4148` at the Suleimanov Fig 2 low-energy longitudinal point.
6. **Thermalization-depth guard.** `magnetic_atmosphere.jl:332-347` correctly uses `τ_eff = ∫√(κ_abs × k_total) dy` (the Sobolev-Mihalas effective optical depth for scattering-dominated atmospheres) as a depth-cutoff in addition to total τ, preventing the LTE bottom BC from being applied above the actual thermalization depth.
7. **Per-mode `B_ν/2` and global `B̄`.** The magnetic Rybicki at `magnetic_atmosphere.jl:482-496, 568-574` correctly splits the unpolarised Planck source per-mode while keeping the global `B̄` summed over both modes — matching the SPW09 + H12 prescription described in the radiative transfer report.
8. **`AtmosphereResult` and `MagneticAtmosphereResult` capture full provenance** of the converged solution (T_eff, g_s, B, θ_B, iteration count, max ΔT/T, grids, T-profile, emergent intensity). Good for downstream analysis.

---

## Closing notes

- **Findings I did not duplicate from reviewer 01:** the `update_atmosphere!` EOS bug (A9 here is a code-smell framing of it), the renderer's silent `θ_B = 0` substitution, and the missing parallelism. Those are architecture-level.
- **Tests** are minimal but well-targeted at the recently-changed magnetic-RT code. Test coverage of the non-magnetic Rybicki step (A5), the `_log_interp` clamp (A24), the `K0_val` fallback (A6), and the `B → 0` limit of `sigma_scat_alpha` (B11) is zero. Reviewer 02 covers this in more depth.
- **NEEDS PAPER VERIFICATION items:** B2 (van Adelsberg & Lai 2006 mode vectors), B7 (PLC2004 Appendix A coefficients), B12 (Ginzburg / P&C 2003 sign of D). All three depend on PDFs that HANDOFF.md flagged as mislabelled or missing.
