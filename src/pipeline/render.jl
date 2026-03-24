#=
End-to-end rendering pipeline (Phase 2 tracer bullet).

Takes NS parameters → renders true-colour and false-colour images.
Uses modified blackbody atmosphere placeholder.
=#

module Renderer

using ..PhysicalConstants: G, c, M_sun, k_B, h, keV
using ..BSkEOS: BSk21_params
using ..TOVSolver: solve_tov
using ..DipoleModel: surface_temperature, surface_Bfield, magnetic_colatitude
using ..BlackbodyAtmosphere: planck_Bnu, emergent_spectrum
using ..SchwarzschildTracer: trace_image, RayResult
using ..CIE_sRGB: load_cie_cmfs, spectrum_to_sRGB, CIE_CMFs

export render_neutron_star, NSParams, RenderResult

"""
Neutron star physical parameters.
"""
struct NSParams
    M_solar::Float64     # mass [M_☉]
    R_km::Float64        # radius [km]
    B_pole::Float64      # polar magnetic field [G]
    T_pole::Float64      # polar temperature [K]
    T_eq::Float64        # equatorial temperature [K]
    obliquity::Float64   # magnetic obliquity [rad]
    inclination::Float64 # observer inclination to spin axis [rad]
    distance_pc::Float64 # distance [pc]
    f_col::Float64       # colour correction factor
end

"""
Result of a rendering.
"""
struct RenderResult
    N::Int
    params::NSParams
    # Per-pixel spectral data: ν_grid and I_ν_obs for each pixel
    ν_grid::Vector{Float64}
    spectra::Array{Vector{Float64}, 2}  # N×N, each element is I_ν at pixel
    # RGB images
    true_colour::Array{NTuple{3,Int}, 2}   # sRGB
    false_colour::Array{NTuple{3,Int}, 2}  # X-ray false colour
end

"""
    render_neutron_star(params, N; n_ν=100) → RenderResult

Full tracer bullet rendering pipeline.
"""
function render_neutron_star(params::NSParams, N::Int;
                             n_ν::Int=100, verbose::Bool=true)
    M = params.M_solar * M_sun
    R = params.R_km * 1e5  # km → cm

    verbose && println("Tracing $(N)×$(N) image...")
    verbose && println("  M = $(params.M_solar) M_☉, R = $(params.R_km) km")
    r_g = 2G * M / c^2
    u = r_g / R
    verbose && println("  u = r_g/R = $(round(u, digits=4)), 1+z = $(round(1/sqrt(1-u), digits=4))")

    # Step 1: Ray trace
    img = trace_image(M, R, params.inclination, N)
    redshift = img.redshift  # 1 + z

    # Step 2: Build frequency grid (0.01 eV to 100 keV, log-spaced)
    # For visible: also need 1.5 - 3.5 eV (350-830 nm)
    E_min = 0.01 * 1.602176634e-12  # 0.01 eV in erg
    E_max = 100.0 * keV             # 100 keV in erg
    ν_grid = [10.0^logν for logν in range(log10(E_min/h), log10(E_max/h), length=n_ν)]

    # Step 3: For each pixel, compute observed spectrum
    spectra = Array{Vector{Float64}}(undef, N, N)
    for j in 1:N, i in 1:N
        ray = img.rays[i, j]
        if !ray.hits
            spectra[i, j] = zeros(n_ν)
            continue
        end

        # Surface properties at hit point
        θ_B = magnetic_colatitude(ray.θ_surface, ray.φ_surface, params.obliquity)
        T_local = surface_temperature(θ_B, params.T_pole, params.T_eq)
        # B_local = surface_Bfield(θ_B, params.B_pole)  # not used in Phase 2

        # Emergent spectrum at surface (in surface rest frame)
        I_ν_em = emergent_spectrum(T_local, ray.cos_α, ν_grid .* redshift;
                                   f_col=params.f_col, p=0.0)

        # Apply redshift: I_ν^obs = (1+z)^{-3} × I_ν'(ν' = (1+z)ν)
        # The factor (1+z)^{-3} comes from I_ν/ν³ being a Lorentz invariant
        spectra[i, j] = I_ν_em ./ redshift^3
    end

    # Step 4: Render true-colour image
    verbose && println("  Computing true-colour rendering...")
    cmfs_path = joinpath(dirname(dirname(@__DIR__)), "refs", "cvrl_cie1931_2deg.csv")
    cmfs = load_cie_cmfs(cmfs_path)

    # Compute average luminance for tone mapping
    Y_total = 0.0
    n_lit = 0
    for j in 1:N, i in 1:N
        if img.rays[i, j].hits && sum(spectra[i, j]) > 0
            # Quick Y estimate from visible-band integral
            Y_total += sum(spectra[i, j])
            n_lit += 1
        end
    end
    L_avg = n_lit > 0 ? Y_total / n_lit : 1.0

    true_colour = Array{NTuple{3,Int}}(undef, N, N)
    for j in 1:N, i in 1:N
        if !img.rays[i, j].hits
            true_colour[i, j] = (0, 0, 0)
        else
            r, g, b = spectrum_to_sRGB(ν_grid, spectra[i, j], cmfs;
                                        exposure=1.0, L_avg=L_avg, key=0.36)
            true_colour[i, j] = (r, g, b)
        end
    end

    # Step 5: False-colour X-ray image
    verbose && println("  Computing false-colour X-ray rendering...")
    false_colour = _false_colour_xray(ν_grid, spectra, N, img)

    verbose && println("  Done.")
    return RenderResult(N, params, ν_grid, spectra, true_colour, false_colour)
end

"""
False-colour X-ray image: map energy bands to RGB.
  Red:   0.1 - 0.5 keV
  Green: 0.5 - 2.0 keV
  Blue:  2.0 - 10.0 keV
"""
function _false_colour_xray(ν_grid, spectra, N, img)
    E_grid = h .* ν_grid ./ keV  # in keV

    fc = Array{NTuple{3,Int}}(undef, N, N)

    # Band boundaries
    r_lo, r_hi = 0.1, 0.5  # keV
    g_lo, g_hi = 0.5, 2.0
    b_lo, b_hi = 2.0, 10.0

    # Compute band fluxes
    R_vals = zeros(N, N)
    G_vals = zeros(N, N)
    B_vals = zeros(N, N)

    for j in 1:N, i in 1:N
        if !img.rays[i, j].hits
            continue
        end
        for k in 1:length(ν_grid)-1
            E = E_grid[k]
            dν = ν_grid[k+1] - ν_grid[k]
            flux = spectra[i, j][k] * dν

            if r_lo <= E < r_hi
                R_vals[i, j] += flux
            elseif g_lo <= E < g_hi
                G_vals[i, j] += flux
            elseif b_lo <= E < b_hi
                B_vals[i, j] += flux
            end
        end
    end

    # Normalize each channel to [0, 255]
    R_max = maximum(R_vals)
    G_max = maximum(G_vals)
    B_max = maximum(B_vals)

    for j in 1:N, i in 1:N
        if !img.rays[i, j].hits
            fc[i, j] = (0, 0, 0)
        else
            r = R_max > 0 ? clamp(round(Int, 255 * R_vals[i, j] / R_max), 0, 255) : 0
            g = G_max > 0 ? clamp(round(Int, 255 * G_vals[i, j] / G_max), 0, 255) : 0
            b = B_max > 0 ? clamp(round(Int, 255 * B_vals[i, j] / B_max), 0, 255) : 0
            fc[i, j] = (r, g, b)
        end
    end

    return fc
end

"""
Save render result as PPM files (portable pixmap — no deps required).
"""
function save_ppm(result::RenderResult, prefix::String)
    N = result.N

    # True colour
    open("$(prefix)_true.ppm", "w") do io
        println(io, "P3")
        println(io, "$N $N")
        println(io, "255")
        for j in N:-1:1  # flip vertically for standard image orientation
            for i in 1:N
                r, g, b = result.true_colour[i, j]
                print(io, "$r $g $b ")
            end
            println(io)
        end
    end

    # False colour
    open("$(prefix)_xray.ppm", "w") do io
        println(io, "P3")
        println(io, "$N $N")
        println(io, "255")
        for j in N:-1:1
            for i in 1:N
                r, g, b = result.false_colour[i, j]
                print(io, "$r $g $b ")
            end
            println(io)
        end
    end
end

end # module
