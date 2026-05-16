#=
Thin newtype unit hierarchy for the NeutronStar user-facing parameter
boundary (bead E10).

Motivation: `NSParams` previously had nine `Float64` fields with six
different unit conventions (R_km in km, T_pole in K, distance_pc in pc,
etc.) — every constructor call was a footgun where `R = 12000` could
plausibly mean meters, km, or cm and the compiler had no way to tell.

This module provides single-field wrappers around `Float64`, each
carrying a fixed CGS-compatible base unit (cm / G / K / g / rad / cm).
Helper constructors `km()`, `gauss()`, `kelvin()`, `M_sun_units()`,
`pc()`, `deg()`, `rad()` etc. perform the conversion at construction.

Design constraint (deliberate, see notes/sessions/2026-05-16-…):
- NO dependency on Unitful.jl — the compile-time cost is large enough
  that the project chose a hand-rolled ~50-line hierarchy instead.
- The newtypes exist ONLY at the user-facing `NSParams` boundary.
  Hot-path solver kwargs (`solve_atmosphere(T_eff::Float64, …)` etc.)
  remain `Float64` in CGS by project convention; callers unwrap via the
  single-field access (`x.cm`, `x.gauss`, …) or via `cgs(x)` for clarity.

The conversion factors `M_sun` and `pc` are pulled from
`PhysicalConstants` to honour the "no magic numbers" CLAUDE.md rule.
=#

module Units

# We deliberately import `PhysicalConstants` as a whole module (not its
# names) because we provide a `pc(x)` helper constructor that would
# collide with `PhysicalConstants.pc` (the numeric constant) if both
# entered this module's namespace by the same short name. The conversion
# factor is therefore referenced fully-qualified below.
import ..PhysicalConstants

export Length, BField, Temperature, Mass, Angle, Distance,
       km, gauss, kelvin, M_sun_units, pc, au, rad, deg, mas, cgs

"Length in centimetres (CGS)."
struct Length
    cm::Float64
end

"Magnetic field strength in Gauss (CGS)."
struct BField
    gauss::Float64
end

"Temperature in Kelvin."
struct Temperature
    kelvin::Float64
end

"Mass in grams (CGS)."
struct Mass
    g::Float64
end

"Angle in radians."
struct Angle
    rad::Float64
end

"Distance in centimetres (CGS). Distinct type from `Length` so that
geometric radii and astronomical distances do not silently substitute
for each other."
struct Distance
    cm::Float64
end

# ── Helper constructors ─────────────────────────────────────────────────
# Each returns the wrapped type with the value converted to the base
# CGS unit. These are plain functions, NOT inner constructors, so users
# always pass through the unit-bearing helper (e.g. `km(12.0)`) rather
# than the raw `Length(12.0)` struct ctor — though the latter is still
# legal for callers that have already done the conversion.

"Kilometres → `Length` (cm)."
km(x::Real)          = Length(Float64(x) * 1.0e5)

"Gauss → `BField`. Identity conversion; provided for symmetry."
gauss(x::Real)       = BField(Float64(x))

"Kelvin → `Temperature`. Identity conversion."
kelvin(x::Real)      = Temperature(Float64(x))

"Solar masses → `Mass` (g). Uses `PhysicalConstants.M_sun`."
M_sun_units(x::Real) = Mass(Float64(x) * PhysicalConstants.M_sun)

"Parsecs → `Distance` (cm). Uses `PhysicalConstants.pc`."
pc(x::Real)          = Distance(Float64(x) * PhysicalConstants.pc)

"Astronomical units → `Distance` (cm). 1 AU = 1.49597870700e13 cm (IAU 2012 exact)."
au(x::Real)          = Distance(Float64(x) * 1.49597870700e13)

"Radians → `Angle`. Identity conversion."
rad(x::Real)         = Angle(Float64(x))

"Degrees → `Angle` (rad)."
deg(x::Real)         = Angle(Float64(x) * π / 180)

"Milli-arcseconds → `Angle` (rad)."
mas(x::Real)         = Angle(Float64(x) * π / (180 * 3.6e6))

# ── Unwrap accessor ─────────────────────────────────────────────────────
# `cgs(x)` is sugar for the single-field accessor. Useful at the
# solver boundary where a `Float64` is required by an internal kwarg.

cgs(x::Length)      = x.cm
cgs(x::BField)      = x.gauss
cgs(x::Temperature) = x.kelvin
cgs(x::Mass)        = x.g
cgs(x::Angle)       = x.rad
cgs(x::Distance)    = x.cm

end # module
