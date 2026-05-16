module NeutronStar

#=
Module dependency order (must match the include sequence below).
Each arrow means "depends on / uses". PhysicalConstants is loaded
first and is a transitive dependency of every other module; arrows
from it are omitted past the first level for readability.

  PhysicalConstants  (loaded first; used by every module below)
        │
        ├──► BSkEOS ───────────► TOVSolver
        │
        ├──► DipoleModel
        │
        ├──► GauntFactor ──┬──► HydrogenOpacity ──┐
        │                  │                      │
        │                  │   MagneticCoulomb ──► MagneticFF
        │                  │                                 │
        │                  │   DielectricTensor ─────────────┤
        │                  │                                 │
        │                  │   BlackbodyAtmosphere ──┬───────┤
        │                  │            │            │       │
        │                  │            │            └──► MagneticModes
        │                  │            │                        │
        │                  └──► AtmosphereStructure ──┐          │
        │                           │      ▲          │          │
        │                           │      └──────────┤          │
        │                           ▼                 │          │
        │                       FeautrierSolver ◄─────┤          │
        │                           │                 │          │
        │                           ▼                 │          │
        │                  TemperatureCorrection ◄────┤          │
        │                           │                 │          │
        │                           ▼                 │          │
        │                       RTAtmosphere ◄────────┤          │
        │                           │                 │          │
        │                           ▼                 │          │
        │                  MagneticAtmosphere ◄───────┴──────────┘
        │                           │
        ├──► SchwarzschildTracer    │
        │           │               │
        ├──► CIE_sRGB               │
        │           │               │
        │           │     AtmosphereGrid ◄── RTAtmosphere, MagneticAtmosphere,
        │           │           │            BlackbodyAtmosphere, FeautrierSolver,
        │           │           │            AtmosphereStructure, GauntFactor
        │           │           ▼
        └───────────┴──────► Renderer ◄── BSkEOS, TOVSolver, DipoleModel,
                                          BlackbodyAtmosphere, SchwarzschildTracer,
                                          CIE_sRGB, AtmosphereGrid

Key invariants:
  * PhysicalConstants must be included first.
  * BlackbodyAtmosphere must come before MagneticModes (uses planck_Bnu)
    and before every atmosphere/* solver.
  * DielectricTensor must come before MagneticModes (polarization weights).
  * MagneticCoulomb must come before MagneticFF (Coulomb log).
  * AtmosphereStructure must come before FeautrierSolver,
    TemperatureCorrection, RTAtmosphere, and MagneticAtmosphere.
  * AtmosphereGrid must come before Renderer (Renderer dispatches on it).

When adding a new module, insert its include() so all modules listed
to its left appear above it. Prefer `using ..ModuleName: symbol` (not
the bare `using ..ModuleName`) so dependencies are explicit at the
import site too — this makes future reorderings safe to verify by
grepping for the imported symbol.
=#

# Physical constants (must be loaded first — all other modules depend on it)
include("constants.jl")
using .PhysicalConstants

# Structured exception types (loaded early so any downstream module may throw them)
include("utils/errors.jl")
using .NSErrors

# Unit newtype hierarchy for the NSParams user-facing boundary (bead E10).
# Loaded right after PhysicalConstants because the helper constructors
# `M_sun_units(x)` and `pc(x)` depend on the numeric constants.
include("utils/units.jl")
using .Units: Length, BField, Temperature, Mass, Angle, Distance,
              km, gauss, kelvin, M_sun_units, pc, au, rad, deg, mas, cgs

# Generic numerical utilities (pure Julia, no internal deps; safe early include)
include("utils/tridiag.jl")
using .Tridiag

# Structured-logging helper for solver entry points (bead E4). Loaded early
# so every downstream solver (rt_atmosphere, magnetic_atmosphere, renderer,
# atmosphere_grid) can wrap its body in `with_solver_logger(verbose) do …`.
include("utils/solver_logging.jl")
using .SolverLogging: with_solver_logger

# Equation of state
include("eos/bsk_eos.jl")
using .BSkEOS

include("eos/tov.jl")
using .TOVSolver

# Surface model
include("surface/dipole.jl")
using .DipoleModel

# Opacities
include("opacity/gaunt_ff.jl")
using .GauntFactor

include("opacity/hydrogen_ff.jl")
using .HydrogenOpacity

include("opacity/coulomb_magnetic.jl")
using .MagneticCoulomb

include("opacity/magnetic_ff.jl")
using .MagneticFF

# Dielectric tensor for normal mode polarisation (must come before magnetic_modes)
include("opacity/dielectric_tensor.jl")
using .DielectricTensor

# Atmosphere (blackbody must come before magnetic_modes which uses planck_Bnu)
include("atmosphere/blackbody.jl")
using .BlackbodyAtmosphere

include("opacity/magnetic_modes.jl")
using .MagneticModes

# Numerical thresholds and damping defaults shared across atmosphere solvers
# (must precede atm_structure / feautrier / rt_atmosphere / magnetic_atmosphere)
include("atmosphere/solver_defaults.jl")
using .SolverDefaults

include("atmosphere/atm_structure.jl")
using .AtmosphereStructure

include("atmosphere/feautrier.jl")
using .FeautrierSolver

include("atmosphere/temp_correction.jl")
using .TemperatureCorrection

include("atmosphere/rt_atmosphere.jl")
using .RTAtmosphere

include("atmosphere/magnetic_atmosphere.jl")
using .MagneticAtmosphere

# GR ray tracing
include("geodesics/schwarzschild.jl")
using .SchwarzschildTracer

# Colorimetry
include("colorimetry/cie_srgb.jl")
using .CIE_sRGB

# Re-export
export PhysicalConstants
export NSErrors, NSError, BadInputError, NumericalError, ConvergenceError, TableOutOfRangeError
# Unit newtypes (bead E10) — user-facing parameter boundary
export Units
export Length, BField, Temperature, Mass, Angle, Distance
export km, gauss, kelvin, M_sun_units, pc, au, rad, deg, mas, cgs
export SolverLogging, with_solver_logger
export BSkEOS, BSk19_params, BSk20_params, BSk21_params
export pressure_of_density, density_of_pressure, energy_density_of_density
export TOVSolver, solve_tov, TOVResult
export DipoleModel, surface_temperature, surface_Bfield, magnetic_colatitude
export GauntFactor, load_gaunt_table, gaunt_ff
export HydrogenOpacity, kappa_ff, sigma_thomson, total_opacity, rosseland_mean
export MagneticFF, cyclotron_freq_e, cyclotron_freq_p
export sigma_ff_alpha, sigma_pp_alpha, sigma_scat_alpha, sigma_total_alpha
export MagneticModes, rosseland_magnetic
export mode_absorption, mode_scattering, mode_opacity, mode_opacity_split, effective_opacity
export SolverDefaults
export TAU_DIFFUSION, TAU_EFFECTIVE_MIN, T_SURFACE_FRAC_MCPHAC,
       Y_MAX_SEMIINFINITE, FLUX_DAMPING_DEFAULT, FLUX_SCALE_CLAMP_LO,
       FLUX_SCALE_CLAMP_HI, DELTA_T_OVER_T_CAP, OPACITY_FLOOR
export AtmosphereStructure, AtmosphereColumn, build_atmosphere, update_atmosphere!
export make_frequency_grid, density_from_PT
export FeautrierSolver, solve_feautrier_all, solve_feautrier_all_adaptive, gauss_legendre_half
export TemperatureCorrection, compute_temperature_correction
export RTAtmosphere, solve_atmosphere, AtmosphereResult, rt_emergent_spectrum
export MagneticAtmosphere, solve_magnetic_atmosphere, MagneticAtmosphereResult, magnetic_emergent_spectrum
export BlackbodyAtmosphere, planck_Bnu, modified_blackbody_intensity, emergent_spectrum
export SchwarzschildTracer, trace_image, RayResult, ImageResult, visible_fraction
export CIE_sRGB, load_cie_cmfs, spectrum_to_XYZ, XYZ_to_linear_sRGB
export linear_sRGB_to_sRGB, tone_map_reinhard, spectrum_to_sRGB

# Pipeline
include("pipeline/atmosphere_grid.jl")
using .AtmosphereGrid
export AtmosphereGrid, AtmosphereSpectrumGrid, AtmosphereGridProvenance, build_atmosphere_grid
export lookup_spectrum, lookup_spectrum_polarized, sum_modes

include("pipeline/render.jl")
using .Renderer
export Renderer, render_neutron_star, render_spectral_cube, NSParams, RenderResult
export SpectralImageCube, render_cube_rgb, save_cube_ppm

end # module NeutronStar
