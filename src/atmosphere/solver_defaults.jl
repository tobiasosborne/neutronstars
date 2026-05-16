"""
    SolverDefaults

Central registry of numerical thresholds and damping parameters used by
the atmosphere solvers. Centralising avoids drift between sites (the
same threshold appearing in 3 files with different values would be a
silent bug).

Edit one value here to retune the solver globally; do not inline these
numbers at use sites.
"""
module SolverDefaults

export TAU_DIFFUSION, TAU_EFFECTIVE_MIN, T_SURFACE_FRAC_MCPHAC,
       Y_MAX_SEMIINFINITE, FLUX_DAMPING_DEFAULT, FLUX_SCALE_CLAMP_LO,
       FLUX_SCALE_CLAMP_HI, DELTA_T_OVER_T_CAP, OPACITY_FLOOR

"Diffusion-depth cutoff for the Feautrier bottom BC (column-integrated optical depth)."
const TAU_DIFFUSION = 80.0

"Thermalization-depth guard: require ∫√(κ_abs·κ_total) dy > this before applying LTE bottom BC."
const TAU_EFFECTIVE_MIN = 10.0

"Empirical non-magnetic surface temperature in units of T_eff (McPHAC; H12 §3.1)."
const T_SURFACE_FRAC_MCPHAC = 0.265

"Semi-infinite-atmosphere column-depth cap [g cm⁻²] (SPW09 scale)."
const Y_MAX_SEMIINFINITE = 1e5

"Default damping factor for the magnetic flux-scaling correction."
const FLUX_DAMPING_DEFAULT = 0.5

"Lower clamp on the per-iteration flux scale (prevents runaway under-correction)."
const FLUX_SCALE_CLAMP_LO = 0.9

"Upper clamp on the per-iteration flux scale (prevents runaway over-correction)."
const FLUX_SCALE_CLAMP_HI = 1.1

"Maximum |ΔT/T| per iteration (corrections rescaled to satisfy this)."
const DELTA_T_OVER_T_CAP = 0.3

"Numerical opacity floor (cm²/g) to avoid divide-by-zero in albedo calculations."
const OPACITY_FLOOR = 1e-30

end # module SolverDefaults
