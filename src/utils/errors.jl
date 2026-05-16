"""
    NSErrors

Structured exception types for the NeutronStar pipeline. Allows callers
to dispatch on failure mode instead of string-matching `ErrorException`s.

Hierarchy:
    NSError
      ├── BadInputError       — user gave physically nonsensical input
      ├── NumericalError      — solver hit a numerical issue (NaN, Inf)
      ├── ConvergenceError    — iterative solver didn't converge within budget
      └── TableOutOfRangeError — interpolation requested outside tabulated range

# Why structured errors?

Phase 4 (Kerr ray tracing on top of TOV) will want to distinguish
recoverable from unrecoverable failures:

  * `BadInputError`        → propagate to user; they retry with better input
  * `NumericalError`       → typically a bug; report with context
  * `ConvergenceError`     → caller may retry with looser tolerance / more iter
  * `TableOutOfRangeError` → caller may fall back to a different model
                             (e.g. fully-ionized atmosphere if partial-ionization
                             table doesn't cover requested (T, ρ))

Without this hierarchy, callers can only `try/catch ErrorException` and
string-match the message — fragile and untyped.
"""
module NSErrors

export NSError, BadInputError, NumericalError, ConvergenceError, TableOutOfRangeError

"""
    NSError <: Exception

Abstract supertype for all NeutronStar pipeline exceptions. Catch this
to handle any structured failure from the package.
"""
abstract type NSError <: Exception end

"""
    BadInputError(msg, [value])

Raised when a user-supplied value is physically nonsensical or outside
the validity range of a model (e.g. central density below 10¹⁰ g/cm³
for a neutron star). The optional `value` field carries the offending
input for diagnostics.
"""
struct BadInputError <: NSError
    msg::String
    value
end
BadInputError(msg::AbstractString) = BadInputError(String(msg), nothing)
Base.showerror(io::IO, e::BadInputError) = print(io, "BadInputError: ", e.msg,
    e.value === nothing ? "" : " (got: $(e.value))")

"""
    NumericalError(msg)

Raised when a solver encounters a NaN, Inf, or other non-finite quantity
mid-computation. Usually indicates a bug or extreme input rather than
something the caller can retry.
"""
struct NumericalError <: NSError
    msg::String
end
Base.showerror(io::IO, e::NumericalError) = print(io, "NumericalError: ", e.msg)

"""
    ConvergenceError(msg, n_iterations, final_residual)

Raised when an iterative solver exhausts its iteration budget without
reaching the requested tolerance. Callers can inspect `n_iterations`
and `final_residual` to decide whether to retry with a looser `rtol`
or a larger iteration cap.
"""
struct ConvergenceError <: NSError
    msg::String
    n_iterations::Int
    final_residual::Float64
end
Base.showerror(io::IO, e::ConvergenceError) = print(io, "ConvergenceError: ", e.msg,
    " (after $(e.n_iterations) iter, final residual $(e.final_residual))")

"""
    TableOutOfRangeError(msg, requested, valid_range)

Raised when an interpolant is queried outside the tabulated domain
(e.g. partial-ionization atmosphere table queried at T or B beyond the
McPHAC grid). `requested` and `valid_range` are intentionally untyped
so callers can pass tuples, scalars, or named ranges as appropriate.
"""
struct TableOutOfRangeError <: NSError
    msg::String
    requested
    valid_range
end
Base.showerror(io::IO, e::TableOutOfRangeError) = print(io, "TableOutOfRangeError: ", e.msg,
    " (requested: $(e.requested), valid: $(e.valid_range))")

end # module NSErrors
