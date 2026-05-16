#=
Shared tridiagonal LU solver used by both the non-magnetic
(TemperatureCorrection) and magnetic (MagneticAtmosphere) Rybicki
temperature corrections.

Two byte-equivalent copies of these helpers previously lived inside
src/atmosphere/temp_correction.jl and src/atmosphere/magnetic_atmosphere.jl.
This module is the canonical version — both consumers `using ..Tridiag`
and call `tridiag_lu_forward!` / `tridiag_lu_back!` from here.

The implementation is a guarded Thomas algorithm: pivots whose absolute
value falls below `PIVOT_TOL` are skipped (treated as zero rows) rather
than dividing by them. This guard exists because Haakonsen Eq. A19-A21
can produce trivially-zero rows at boundaries that should propagate as
no-ops instead of NaNs.
=#

module Tridiag

export tridiag_lu_forward!, tridiag_lu_back!

"Pivot magnitude below this is treated as zero (skip the row)."
const PIVOT_TOL = 1e-30

"""
    tridiag_lu_forward!(sub, diag, sup, d) -> (lu_diag, lu_sup, lu_rhs)

LU forward sweep for a scalar tridiagonal system with sub-diagonal `sub`
(length N-1), diagonal `diag` (length N), super-diagonal `sup` (length
N-1) and right-hand side `d` (length N).

The inputs are not mutated; fresh `lu_diag`, `lu_sup`, `lu_rhs` arrays
are returned (matching the contract of the original
`_tridiag_lu_forward` helpers in temp_correction.jl /
magnetic_atmosphere.jl).

Rows whose pivot magnitude is below `PIVOT_TOL` are skipped, which lets
boundary rows with structurally-zero diagonal pass through cleanly.
"""
function tridiag_lu_forward!(sub::AbstractVector, diag::AbstractVector,
                              sup::AbstractVector, d::AbstractVector)
    N = length(diag)
    lu_diag = copy(diag)
    lu_sup = copy(sup)
    lu_rhs = copy(d)

    for i in 2:N
        if abs(lu_diag[i-1]) < PIVOT_TOL
            continue
        end
        m = sub[i-1] / lu_diag[i-1]
        lu_diag[i] -= m * lu_sup[i-1]
        lu_rhs[i] -= m * lu_rhs[i-1]
    end

    return lu_diag, lu_sup, lu_rhs
end

"""
    tridiag_lu_back!(lu_diag, lu_sup, lu_rhs) -> x

Back-substitution on the factored arrays returned by
`tridiag_lu_forward!`. Returns a freshly-allocated solution vector `x`
of length `length(lu_diag)`.

Rows whose pivot magnitude is below `PIVOT_TOL` leave the corresponding
`x[i]` at zero, matching the original helpers' behaviour.

The output type promotes the input element types so this composes through
ForwardDiff `Dual` numbers (the AD-percolation follow-up to bead E6).
"""
function tridiag_lu_back!(lu_diag::AbstractVector, lu_sup::AbstractVector,
                           lu_rhs::AbstractVector)
    N = length(lu_diag)
    R = promote_type(eltype(lu_diag), eltype(lu_sup), eltype(lu_rhs))
    x = zeros(R, N)
    if abs(lu_diag[N]) > PIVOT_TOL
        x[N] = lu_rhs[N] / lu_diag[N]
    end
    for i in N-1:-1:1
        if abs(lu_diag[i]) > PIVOT_TOL
            x[i] = (lu_rhs[i] - lu_sup[i] * x[i+1]) / lu_diag[i]
        end
    end
    return x
end

end # module
