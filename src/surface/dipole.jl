#=
Neutron star surface model: dipole magnetic field and Greenstein-Hartke T(θ).

Source: Greenstein & Hartke, ApJ 271, 283 (1983).
Local: refs/greenstein_hartke_1983_surface_T.pdf
=#

module DipoleModel

export surface_temperature, surface_Bfield, magnetic_colatitude

"""
    magnetic_colatitude(θ_geo, φ_geo, obliquity) → θ_B

Magnetic colatitude at geographic position (θ_geo, φ_geo) for a dipole
tilted by `obliquity` from the spin axis.

The magnetic axis is in the (sin α, 0, cos α) direction where α = obliquity.
"""
function magnetic_colatitude(θ_geo::Float64, φ_geo::Float64,
                              obliquity::Float64)::Float64
    # Dot product of surface normal with magnetic axis
    cos_θB = cos(θ_geo) * cos(obliquity) +
             sin(θ_geo) * cos(φ_geo) * sin(obliquity)
    return acos(clamp(cos_θB, -1.0, 1.0))
end

"""
    surface_temperature(θ_B, T_pole, T_eq) → T [K]

Greenstein-Hartke temperature model.
T(θ_B) = T_eq + (T_pole - T_eq) cos²(θ_B)

Source: Greenstein & Hartke (1983) Eq. 1.
"""
function surface_temperature(θ_B::Float64, T_pole::Float64,
                              T_eq::Float64)::Float64
    @assert T_pole > 0 "T_pole must be positive"
    @assert T_eq > 0 "T_eq must be positive"
    @assert T_pole >= T_eq "T_pole should be ≥ T_eq"

    T = T_eq + (T_pole - T_eq) * cos(θ_B)^2
    @assert T > 0 "Temperature must be positive"
    return T
end

"""
    surface_Bfield(θ_B, B_pole) → |B| [Gauss]

Magnetic field magnitude at magnetic colatitude θ_B for a centred dipole.
|B| = (B_pole/2) √(1 + 3cos²θ_B)

Source: Standard dipole field. Greenstein & Hartke (1983).
"""
function surface_Bfield(θ_B::Float64, B_pole::Float64)::Float64
    @assert B_pole >= 0 "B_pole must be non-negative"
    return 0.5 * B_pole * sqrt(1.0 + 3.0 * cos(θ_B)^2)
end

end # module
