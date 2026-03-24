module NeutronStar

# Physical constants (must be loaded first — all other modules depend on it)
include("constants.jl")
using .PhysicalConstants

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

include("atmosphere/atm_structure.jl")
using .AtmosphereStructure

include("atmosphere/feautrier.jl")
using .FeautrierSolver

include("atmosphere/temp_correction.jl")
using .TemperatureCorrection

include("atmosphere/rt_atmosphere.jl")
using .RTAtmosphere

# GR ray tracing
include("geodesics/schwarzschild.jl")
using .SchwarzschildTracer

# Colorimetry
include("colorimetry/cie_srgb.jl")
using .CIE_sRGB

# Re-export
export PhysicalConstants
export BSkEOS, BSk19_params, BSk20_params, BSk21_params
export pressure_of_density, density_of_pressure, energy_density_of_density
export TOVSolver, solve_tov, TOVResult
export DipoleModel, surface_temperature, surface_Bfield, magnetic_colatitude
export GauntFactor, load_gaunt_table, gaunt_ff
export HydrogenOpacity, kappa_ff, sigma_thomson, total_opacity, rosseland_mean
export MagneticFF, cyclotron_freq_e, cyclotron_freq_p, sigma_ff_alpha, sigma_total_alpha
export MagneticModes, rosseland_magnetic, mode_opacity, effective_opacity
export AtmosphereStructure, AtmosphereColumn, build_atmosphere, update_atmosphere!
export make_frequency_grid
export FeautrierSolver, solve_feautrier_all, gauss_legendre_half
export TemperatureCorrection, compute_temperature_correction
export RTAtmosphere, solve_atmosphere, AtmosphereResult, rt_emergent_spectrum
export BlackbodyAtmosphere, planck_Bnu, modified_blackbody_intensity, emergent_spectrum
export SchwarzschildTracer, trace_image, RayResult, ImageResult, visible_fraction
export CIE_sRGB, load_cie_cmfs, spectrum_to_XYZ, XYZ_to_linear_sRGB
export linear_sRGB_to_sRGB, tone_map_reinhard, spectrum_to_sRGB

# Pipeline
include("pipeline/render.jl")
using .Renderer
export Renderer, render_neutron_star, NSParams, RenderResult

end # module NeutronStar
