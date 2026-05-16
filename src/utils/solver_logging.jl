#=
Solver logging helper (bead E4).

Centralises the `verbose::Bool` → logger contract used by every long-
running solver in the project. Previously each entry point sprinkled
`verbose && @printf(...)` / `verbose && println(...)` calls throughout
its body, which:

  * put the formatting strings and the on/off behaviour in different
    places (and required the caller to thread `verbose` through every
    helper that wanted to log),
  * gave downstream callers no way to redirect output (log file, JSON
    capture, test instrumentation) without monkey-patching `Base.stdout`.

This module replaces that pattern with Julia's standard `Logging`
macros (`@info`, `@debug`, `@warn`). Each solver entry point wraps its
body in `with_solver_logger(verbose) do … end`, which routes to
`NullLogger` when `verbose=false` and to the active logger otherwise.
Per-iteration progress lines are demoted to `@debug` so they stay
silent at the default `Info` level and light up when the user sets
`JULIA_DEBUG=NeutronStar`.
=#

module SolverLogging

using Logging

export with_solver_logger

"""
    with_solver_logger(f, verbose::Bool)

Run `f()` under either the active logger (when `verbose=true`) or
`NullLogger` (when `verbose=false`). The function comes first so the
`do`-block form works:

    with_solver_logger(verbose) do
        @info "starting solve"
        ...
    end

(Julia's `do` syntax expands `g(args...) do x; ...; end` into
`g(x -> ..., args...)`, i.e. it injects the closure as the *first*
positional argument — so the callable must come before `verbose` in
the signature.)

The wrapper preserves the existing `verbose::Bool` kwarg contract:
`verbose=false` is fully silent; `verbose=true` emits Info+Warn to
the default logger (which downstream callers can swap via
`with_logger`). Per-iteration progress should use `@debug` so it stays
silent under default Info-level filtering and lights up with
`JULIA_DEBUG=NeutronStar`.
"""
function with_solver_logger(f, verbose::Bool)
    return verbose ? f() : with_logger(f, NullLogger())
end

end # module
