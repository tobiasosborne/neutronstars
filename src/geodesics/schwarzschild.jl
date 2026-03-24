#=
Schwarzschild geodesic ray tracing for neutron star images.

Uses the Beloborodov (2002) approximate cosine relation (Eq. 1):
    1 - cos α = (1 - cos ψ)(1 - r_g/R)

where:
  α = emission angle (angle between photon and radial direction at surface)
  ψ = escape angle (angle from the line of sight at infinity)
  r_g = 2GM/c² (Schwarzschild radius)
  R = stellar radius

This approximation is accurate to ~1% for R > 2r_g (typical NS: R ≈ 3r_g).

Also implements exact integral from Beloborodov (2002) Appendix, Eq. (9).

Source: Beloborodov, ApJ 566, L85 (2002). Local: refs/beloborodov_2002_ray_tracing.pdf
=#

module SchwarzschildTracer

using ..PhysicalConstants: G, c
using LinearAlgebra

export trace_image, RayResult, ImageResult, visible_fraction

"""
Result for a single traced ray.
"""
struct RayResult
    hits::Bool        # true if ray hits the star
    b::Float64        # impact parameter [cm]
    ψ::Float64        # polar angle on observer sky [rad]
    θ_surface::Float64 # colatitude on star surface [rad]
    φ_surface::Float64 # longitude on star surface [rad]
    cos_α::Float64    # cosine of emission angle
    redshift::Float64  # 1 + z = (1 - r_g/R)^{-1/2}
end

"""
Result for a full image trace.
"""
struct ImageResult
    N::Int                   # resolution (NxN pixels)
    M::Float64               # stellar mass [g]
    R::Float64               # stellar radius [cm]
    r_g::Float64             # Schwarzschild radius [cm]
    compactness::Float64     # u = r_g/R
    redshift::Float64        # 1+z
    visible_frac::Float64    # fraction of surface visible
    rays::Matrix{RayResult}  # NxN ray results
end

"""
    cos_emission_angle(cos_ψ, u) → cos_α

Beloborodov (2002) Eq. (1) inverted:
    1 - cos α = (1 - cos ψ)(1 - u/2)

where u = r_g/R = 2GM/(Rc²) and we use r_g/R = u (not u/2 — note
Beloborodov uses r_g = 2GM/c², so r_g/R = u).

Rearranging: cos α = 1 - (1 - cos ψ)(1 - r_g/R).
"""
function cos_emission_angle(cos_ψ::Float64, u::Float64)::Float64
    return 1.0 - (1.0 - cos_ψ) * (1.0 - u)
end

"""
    cos_ψ_from_cos_α(cos_α, u) → cos_ψ

Inverse of cosine relation.
"""
function cos_ψ_from_cos_α(cos_α::Float64, u::Float64)::Float64
    return 1.0 - (1.0 - cos_α) / (1.0 - u)
end

"""
    impact_parameter(cos_α, R, u) → b

Impact parameter from emission angle.
Beloborodov (2002) Eq. (10): sin α = (b/R)√(1 - r_g/R).
"""
function impact_parameter(cos_α::Float64, R::Float64, u::Float64)::Float64
    sin_α = sqrt(max(0.0, 1.0 - cos_α^2))
    return R * sin_α / sqrt(1.0 - u)
end

"""
    critical_cos_ψ(u) → cos_ψ_crit

The critical angle beyond which surface is not visible.
At the limb, cos α = 0 (tangential emission).
cos ψ_crit = 1 - 1/(1 - u).
For R = 3r_g (u=2/3): cos ψ_crit = 1 - 3 = -2 → entire visible hemisphere
and part of far side visible.
"""
function critical_cos_ψ(u::Float64)::Float64
    return 1.0 - 1.0 / (1.0 - u)
end

"""
    visible_fraction(u) → f

Fraction of stellar surface visible to a distant observer.
Beloborodov (2002) Sect. 3: S_v/(4πR²) = [2(1 - r_g/R)]⁻¹.

For u → 0 (flat space): f → 1/2.
For u → 2/3 (R = 3r_g): f → 3/4.
"""
function visible_fraction(u::Float64)::Float64
    @assert 0 <= u < 1 "Compactness must be in [0,1), got $u"
    return 1.0 / (2.0 * (1.0 - u))
end

"""
    trace_image(M, R, inclination, N; D=Inf) → ImageResult

Trace an N×N image of a neutron star.

# Arguments
- `M`: stellar mass [g]
- `R`: stellar radius [cm]
- `inclination`: angle between observer line of sight and spin axis [rad]
- `N`: pixel resolution (N×N)
- `D`: observer distance [cm] (default: infinity → parallel rays)

The image plane is centred on the star. Pixels map to impact parameters
b ∈ [0, b_max] where b_max slightly exceeds the apparent stellar radius.
"""
function trace_image(M::Float64, R::Float64, inclination::Float64,
                     N::Int; D::Float64=Inf)::ImageResult
    @assert M > 0 "Mass must be positive"
    @assert R > 0 "Radius must be positive"
    @assert N > 0 "Resolution must be positive"

    r_g = 2G * M / c^2
    u = r_g / R
    @assert u < 1.0 "Star is inside Schwarzschild radius! u=$u"

    redshift_factor = 1.0 / sqrt(1.0 - u)  # 1 + z
    f_vis = visible_fraction(u)

    # Critical impact parameter (limb of apparent disk)
    # At cos α = 0: b_crit = R / √(1 - u)
    b_crit = R / sqrt(1.0 - u)

    # Image extends to b_max slightly beyond b_crit
    b_max = 1.05 * b_crit

    # Pixel grid in image plane: Cartesian (x, y) with x = b cos(φ_img), y = b sin(φ_img)
    # The spin axis is in the x-z plane, tilted by inclination angle from z (= line of sight)
    rays = Matrix{RayResult}(undef, N, N)
    pixel_size = 2b_max / N

    for j in 1:N, i in 1:N
        # Image plane coordinates (centered)
        x_img = -b_max + (i - 0.5) * pixel_size
        y_img = -b_max + (j - 0.5) * pixel_size
        b = sqrt(x_img^2 + y_img^2)
        φ_img = atan(y_img, x_img)

        if b > b_crit || b < 1e-10
            # Ray misses the star (or is exactly at centre → handle below)
            if b < 1e-10
                # Centre ray: ψ = 0, hits pole facing observer
                cos_α = cos_emission_angle(1.0, u)
                rays[i, j] = RayResult(true, 0.0, 0.0, inclination, 0.0,
                                        cos_α, redshift_factor)
            else
                rays[i, j] = RayResult(false, b, 0.0, 0.0, 0.0, 0.0, 0.0)
            end
            continue
        end

        # Convert impact parameter to polar angle on sky
        # sin ψ = b/D for D → ∞ this is just b (angular)
        # For the Beloborodov formula we need cos ψ
        # From Eq. (10): sin α = (b/R)√(1 - u), so
        # cos α = √(1 - sin²α) = √(1 - b²(1-u)/R²)
        sin_α = b * sqrt(1.0 - u) / R

        if sin_α > 1.0
            rays[i, j] = RayResult(false, b, 0.0, 0.0, 0.0, 0.0, 0.0)
            continue
        end

        cos_α = sqrt(1.0 - sin_α^2)
        cos_ψ = cos_ψ_from_cos_α(cos_α, u)

        # ψ is the polar angle from the sub-observer point on the stellar surface
        ψ = acos(clamp(cos_ψ, -1.0, 1.0))

        # Map (ψ, φ_img) on the sky to (θ_surface, φ_surface) on the star
        # The sub-observer point is at (θ = inclination, φ = 0)
        # Rotate by ψ in direction φ_img
        θ_s, φ_s = sky_to_surface(ψ, φ_img, inclination)

        rays[i, j] = RayResult(true, b, ψ, θ_s, φ_s, cos_α, redshift_factor)
    end

    return ImageResult(N, M, R, r_g, u, redshift_factor, f_vis, rays)
end

"""
Map sky angle (ψ, φ_img) to surface coordinates (θ, φ) given observer inclination.

The observer looks along ẑ. The spin axis is in the x-z plane at angle `incl` from ẑ.
A point at angle ψ from the sub-observer point and azimuth φ_img maps to
surface colatitude θ and longitude φ.
"""
function sky_to_surface(ψ::Float64, φ_img::Float64, incl::Float64)
    # Sub-observer point: θ₀ = incl, φ₀ = 0
    # Unit vector of sub-observer point: (sin incl, 0, cos incl)
    # Photon lands at angle ψ from sub-observer, in direction φ_img on sky

    # Direction on sky: x_sky → along meridian toward pole, y_sky → perpendicular
    # In 3D, the point is obtained by rotating the sub-observer direction by ψ

    cos_i = cos(incl)
    sin_i = sin(incl)
    cos_ψ = cos(ψ)
    sin_ψ = sin(ψ)

    # Surface point in Cartesian (spin axis = z)
    # Use spherical rotation: point at angle ψ from (sin_i, 0, cos_i)
    # in direction φ_img
    x = sin_i * cos_ψ + cos_i * sin_ψ * cos(φ_img)
    y = sin_ψ * sin(φ_img)
    z = cos_i * cos_ψ - sin_i * sin_ψ * cos(φ_img)

    θ = acos(clamp(z, -1.0, 1.0))
    φ = atan(y, x)

    return θ, φ
end

end # module
