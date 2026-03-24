#=
Modified blackbody atmosphere placeholder (Phase 2 tracer bullet).

I_ν = f_col⁴ × B_ν(f_col × T, ν) × cos(θ_e)^p

where f_col is the colour correction factor and p is a limb-darkening exponent.

APPROXIMATION: This is a placeholder for self-consistent radiative transfer.
See af node 1.3.1.

Source: Suleimanov+ (2009) for f_col values.
Planck function: standard, B_ν = 2hν³/c² × 1/(exp(hν/kT) - 1).
=#

module BlackbodyAtmosphere

using ..PhysicalConstants: h, c, k_B

export planck_Bnu, modified_blackbody_intensity, emergent_spectrum

"""
    planck_Bnu(ν, T) → B_ν [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹]

Planck function. Source: standard thermodynamics.
"""
function planck_Bnu(ν::Float64, T::Float64)::Float64
    @assert ν > 0 "Frequency must be positive"
    @assert T > 0 "Temperature must be positive"

    x = h * ν / (k_B * T)

    if x > 500.0
        return 0.0  # exp overflow protection
    end

    return 2h * ν^3 / c^2 / (exp(x) - 1.0)
end

"""
    modified_blackbody_intensity(ν, T, cos_θe; f_col=1.0, p=0.0) → I_ν

Modified blackbody with colour correction and limb darkening.

# Arguments
- `ν`: photon frequency [Hz]
- `T`: local effective temperature [K]
- `cos_θe`: cosine of emission angle
- `f_col`: colour correction factor (default 1.0 = pure blackbody)
- `p`: limb-darkening exponent (default 0.0 = isotropic)
"""
function modified_blackbody_intensity(ν::Float64, T::Float64,
                                       cos_θe::Float64;
                                       f_col::Float64=1.0,
                                       p::Float64=0.0)::Float64
    @assert cos_θe >= 0 "cos_θe must be ≥ 0 (visible hemisphere)"
    @assert f_col > 0 "f_col must be positive"

    B = planck_Bnu(ν, f_col * T)
    limb = cos_θe > 0 ? cos_θe^p : 0.0

    return f_col^4 * B * limb
end

"""
    emergent_spectrum(T, cos_θe, ν_grid; f_col=1.0, p=0.0) → I_ν[]

Compute spectrum I_ν at each frequency in ν_grid.
"""
function emergent_spectrum(T::Float64, cos_θe::Float64,
                           ν_grid::AbstractVector{Float64};
                           f_col::Float64=1.0,
                           p::Float64=0.0)::Vector{Float64}
    return [modified_blackbody_intensity(ν, T, cos_θe; f_col=f_col, p=p)
            for ν in ν_grid]
end

end # module
