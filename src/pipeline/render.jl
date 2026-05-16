#=
End-to-end rendering pipeline (Phase 2 tracer bullet).

Takes NS parameters → renders true-colour and false-colour images.
Uses modified blackbody atmosphere placeholder.
=#

module Renderer

using Printf
using Logging
using ..PhysicalConstants: G, c, M_sun, k_B, h, keV, pc
using ..Units: Length, BField, Temperature, Mass, Angle, Distance
using ..BSkEOS: BSk21_params
using ..TOVSolver: solve_tov
using ..DipoleModel: surface_temperature, surface_Bfield, magnetic_colatitude
using ..BlackbodyAtmosphere: planck_Bnu, emergent_spectrum
using ..SchwarzschildTracer: trace_image, RayResult
using ..CIE_sRGB: load_cie_cmfs, spectrum_to_sRGB, CIE_CMFs
using ..AtmosphereGrid: AtmosphereSpectrumGrid, lookup_spectrum
using ..SolverLogging: with_solver_logger

export render_neutron_star, render_spectral_cube, NSParams, RenderResult
export SpectralImageCube, render_cube_rgb, save_cube_ppm

"""
Neutron star physical parameters.

All fields except `f_col` are dimensional newtypes (see `Units` module,
bead E10) so the compiler refuses constructor calls that mix units.
Build them with the unit-bearing helpers:

    params = NSParams(
        M_sun_units(1.4),   # mass [M_☉]      → Mass        (grams CGS)
        km(12.0),           # radius [km]     → Length      (cm CGS)
        gauss(1e12),        # polar B [G]     → BField
        kelvin(1.5e6),      # T_pole [K]      → Temperature
        kelvin(3e5),        # T_eq [K]        → Temperature
        rad(0.3),           # obliquity [rad] → Angle
        rad(π/3),           # inclination     → Angle
        pc(100.0),          # distance [pc]   → Distance    (cm CGS)
        1.0,                # f_col (dimensionless Float64)
    )

Fields:
- `M`           : mass (`Mass`; internal unit g). User builds with `M_sun_units`.
- `R`           : radius (`Length`; internal unit cm). User builds with `km`.
- `B_pole`      : polar magnetic field (`BField`; G).
- `T_pole`      : polar effective temperature (`Temperature`; K).
- `T_eq`        : equatorial effective temperature (`Temperature`; K).
- `obliquity`   : magnetic obliquity (`Angle`; rad).
- `inclination` : observer inclination to spin axis (`Angle`; rad).
- `distance`    : distance (`Distance`; internal unit cm). User builds with `pc` or `au`.
- `f_col`       : LEGACY colour-correction factor for the modified-blackbody
                  atmosphere placeholder. **Only used by `render_neutron_star`.**
                  `render_spectral_cube` ignores this field — colour correction
                  is already baked into the pre-computed atmosphere grid spectra,
                  so applying `f_col` again would be double-counting. A non-unity
                  value combined with `render_spectral_cube` triggers a `@warn`.

The legacy modified-blackbody spectral-hardening parameter `p` from
`BlackbodyAtmosphere.emergent_spectrum` is hardcoded to `0.0` in
`render_neutron_star` and is not exposed here.
"""
struct NSParams
    M::Mass              # mass; .g gives grams (CGS)
    R::Length            # radius; .cm gives cm (CGS)
    B_pole::BField       # polar magnetic field; .gauss gives G
    T_pole::Temperature  # polar temperature; .kelvin gives K
    T_eq::Temperature    # equatorial temperature; .kelvin gives K
    obliquity::Angle     # magnetic obliquity; .rad gives rad
    inclination::Angle   # observer inclination to spin axis; .rad gives rad
    distance::Distance   # distance; .cm gives cm (CGS)
    f_col::Float64       # colour correction factor (legacy; render_neutron_star only)
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
    M = params.M.g          # grams (already converted at construction)
    R = params.R.cm         # cm    (already converted at construction)

    return with_solver_logger(verbose) do
    @info "Tracing $(N)×$(N) image..."
    @info "  M = $(M / M_sun) M_☉, R = $(R / 1e5) km"
    r_g = 2G * M / c^2
    u = r_g / R
    @info "  u = r_g/R = $(round(u, digits=4)), 1+z = $(round(1/sqrt(1-u), digits=4))"

    # Step 1: Ray trace
    img = trace_image(M, R, params.inclination.rad, N)
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
        θ_B = magnetic_colatitude(ray.θ_surface, ray.φ_surface, params.obliquity.rad)
        T_local = surface_temperature(θ_B, params.T_pole.kelvin, params.T_eq.kelvin)
        # B_local = surface_Bfield(θ_B, params.B_pole.gauss)  # not used in Phase 2

        # Emergent spectrum at surface (in surface rest frame)
        I_ν_em = emergent_spectrum(T_local, ray.cos_α, ν_grid .* redshift;
                                   f_col=params.f_col, p=0.0)

        # Apply redshift: I_ν^obs = (1+z)^{-3} × I_ν'(ν' = (1+z)ν)
        # The factor (1+z)^{-3} comes from I_ν/ν³ being a Lorentz invariant
        spectra[i, j] = I_ν_em ./ redshift^3
    end

    # Step 4: Render true-colour image
    @info "  Computing true-colour rendering..."
    cmfs_path = joinpath(@__DIR__, "..", "data", "cvrl_cie1931_2deg.csv")
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
    @info "  Computing false-colour X-ray rendering..."
    false_colour = _false_colour_xray(ν_grid, spectra, N, img)

    @info "  Done."
    return RenderResult(N, params, ν_grid, spectra, true_colour, false_colour)
    end  # with_solver_logger
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

"""
Spectral image cube: the primary data product.
I(x, y, ν) — specific intensity at each pixel and frequency.
"""
struct SpectralImageCube
    nx::Int; ny::Int; nν::Int
    ν_grid::Vector{Float64}           # observer-frame frequencies [Hz]
    pixel_scale::Float64              # angular size per pixel [rad]
    I::Array{Float64, 3}             # (nx, ny, nν) specific intensity [erg/s/cm²/Hz/sr]
    hit_fraction::Matrix{Float64}     # fraction of pixel hitting surface
    mean_redshift::Matrix{Float64}
    mean_T::Matrix{Float64}
    mean_B::Matrix{Float64}
    params::NSParams
    atmosphere_id::String
end

"""
    render_spectral_cube(params, atm_grid, N; verbose=true) → SpectralImageCube

Render a neutron star using pre-computed atmosphere spectra.
Replaces the modified-blackbody placeholder with real RT solutions.

For each pixel:
  1. Ray trace → surface coordinates, emission angle, redshift
  2. Surface model → local T_eff, local B
  3. Look up atmosphere spectrum I_ν(cos θ_e) from pre-computed grid
  4. Apply relativistic transform: I_ν^obs = I_ν'(ν') / (1+z)³
"""
function render_spectral_cube(params::NSParams,
                               atm_grid::AtmosphereSpectrumGrid,
                               N::Int; verbose::Bool=true)
    # Guard against the post-D1 footgun: NSParams.f_col is a legacy knob for
    # the modified-blackbody path (render_neutron_star). The atmosphere-grid
    # path bakes colour correction into the pre-computed spectra, so f_col
    # here would be silently ignored. This @warn fires regardless of `verbose`
    # because it indicates a real bug in the caller's setup.
    if params.f_col != 1.0
        @warn "NSParams.f_col is ignored by render_spectral_cube (only the legacy render_neutron_star uses it); colour correction is already encoded in the atmosphere grid spectra" f_col=params.f_col
    end

    M_cgs = params.M.g       # grams (already converted at construction)
    R = params.R.cm          # cm    (already converted at construction)

    return with_solver_logger(verbose) do
    @info "Rendering $(N)×$(N) spectral cube with real atmosphere spectra..."
    r_g = 2G * M_cgs / c^2
    u = r_g / R
    z_surf = 1.0 / sqrt(1.0 - u) - 1.0
    redshift = 1.0 + z_surf
    @info @sprintf("  u=r_g/R=%.4f, 1+z=%.4f", u, redshift)

    # Ray trace
    @info "  Ray tracing..."
    img = trace_image(M_cgs, R, params.inclination.rad, N)

    # Use the atmosphere grid's frequency grid (observer frame)
    # ν_obs = ν_em / (1+z), so ν_em = ν_obs × (1+z)
    ν_grid_em = atm_grid.ν_grid
    nν = length(ν_grid_em)
    ν_grid_obs = ν_grid_em ./ redshift

    # Allocate cube
    I_cube = zeros(N, N, nν)
    hit_frac = zeros(N, N)
    mean_z = zeros(N, N)
    mean_T = zeros(N, N)
    mean_B_arr = zeros(N, N)

    @info "  Computing spectra for $(N*N) pixels, $(nν) frequencies..."
    n_hit = 0
    for j in 1:N, i in 1:N
        ray = img.rays[i, j]
        if !ray.hits
            continue
        end
        n_hit += 1

        # Surface properties at hit point
        θ_B = magnetic_colatitude(ray.θ_surface, ray.φ_surface, params.obliquity.rad)
        T_local = surface_temperature(θ_B, params.T_pole.kelvin, params.T_eq.kelvin)
        B_local = surface_Bfield(θ_B, params.B_pole.gauss)

        hit_frac[i, j] = 1.0
        mean_z[i, j] = redshift
        mean_T[i, j] = T_local
        mean_B_arr[i, j] = B_local

        # Emission angle: cos θ_e from ray tracer
        cos_θe = ray.cos_α

        # Look up atmosphere spectrum at emission frequencies.
        # IMPORTANT (post-D1): θ_B is the angle between B and the surface
        # normal at the hit point and is the dominant variable for the
        # X/O-mode opacity ratio. Pre-D1 this was hardcoded to 0 inside
        # lookup_spectrum, silently using pole-on opacities for every pixel.
        I_em = lookup_spectrum(atm_grid, T_local, B_local, θ_B,
                                ν_grid_em, cos_θe)

        # Relativistic transform: I_ν^obs = I_ν'(ν') / (1+z)³
        for k in 1:nν
            I_cube[i, j, k] = I_em[k] / redshift^3
        end
    end

    @info @sprintf("  %d/%d pixels hit surface (%.0f%%)",
                   n_hit, N*N, 100.0 * n_hit / (N*N))

    # Pixel angular scale (distance is already in cm; ratio is dimensionless)
    pixel_scale = 2.0 * R / params.distance.cm / N  # approximate

    @info "  Spectral cube complete."
    return SpectralImageCube(N, N, nν, ν_grid_obs, pixel_scale,
                              I_cube, hit_frac, mean_z, mean_T, mean_B_arr,
                              params, "AtmosphereGrid_v1")
    end  # with_solver_logger
end

"""
    render_cube_rgb(cube; cmfs_path) → (true_colour, false_colour)

Post-process a SpectralImageCube into RGB images.
"""
function render_cube_rgb(cube::SpectralImageCube;
                          cmfs_path::String="")
    N = cube.nx
    ν_grid = cube.ν_grid

    if cmfs_path == ""
        cmfs_path = joinpath(@__DIR__, "..", "data", "cvrl_cie1931_2deg.csv")
    end
    cmfs = load_cie_cmfs(cmfs_path)

    # Compute average luminance
    Y_total = 0.0; n_lit = 0
    for j in 1:N, i in 1:N
        if cube.hit_fraction[i, j] > 0 && sum(cube.I[i, j, :]) > 0
            Y_total += sum(cube.I[i, j, :]); n_lit += 1
        end
    end
    L_avg = n_lit > 0 ? Y_total / n_lit : 1.0

    # True colour via CIE 1931
    true_colour = Array{NTuple{3,Int}}(undef, N, N)
    for j in 1:N, i in 1:N
        if cube.hit_fraction[i, j] == 0
            true_colour[i, j] = (0, 0, 0)
        else
            r, g, b = spectrum_to_sRGB(ν_grid, cube.I[i, j, :], cmfs;
                                        exposure=1.0, L_avg=L_avg, key=0.36)
            true_colour[i, j] = (r, g, b)
        end
    end

    # False colour X-ray
    false_colour = _false_colour_xray(ν_grid, cube, N)

    return true_colour, false_colour
end

"""False-colour X-ray from SpectralImageCube."""
function _false_colour_xray(ν_grid, cube::SpectralImageCube, N)
    E_grid = h .* ν_grid ./ keV
    fc = Array{NTuple{3,Int}}(undef, N, N)
    R_vals = zeros(N, N); G_vals = zeros(N, N); B_vals = zeros(N, N)

    for j in 1:N, i in 1:N
        cube.hit_fraction[i, j] == 0 && continue
        for k in 1:length(ν_grid)-1
            E = E_grid[k]; dν = ν_grid[k+1] - ν_grid[k]
            flux = cube.I[i, j, k] * dν
            if 0.1 <= E < 0.5;     R_vals[i, j] += flux
            elseif 0.5 <= E < 2.0; G_vals[i, j] += flux
            elseif 2.0 <= E < 10.0; B_vals[i, j] += flux
            end
        end
    end

    Rm = maximum(R_vals); Gm = maximum(G_vals); Bm = maximum(B_vals)
    for j in 1:N, i in 1:N
        if cube.hit_fraction[i, j] == 0
            fc[i, j] = (0, 0, 0)
        else
            r = Rm > 0 ? clamp(round(Int, 255 * R_vals[i, j] / Rm), 0, 255) : 0
            g = Gm > 0 ? clamp(round(Int, 255 * G_vals[i, j] / Gm), 0, 255) : 0
            b = Bm > 0 ? clamp(round(Int, 255 * B_vals[i, j] / Bm), 0, 255) : 0
            fc[i, j] = (r, g, b)
        end
    end
    return fc
end

"""Save SpectralImageCube renderings as PPM."""
function save_cube_ppm(true_colour, false_colour, N::Int, prefix::String)
    for (img, suffix) in [(true_colour, "true"), (false_colour, "xray")]
        open("$(prefix)_$(suffix).ppm", "w") do io
            println(io, "P3"); println(io, "$N $N"); println(io, "255")
            for j in N:-1:1
                for i in 1:N
                    r, g, b = img[i, j]
                    print(io, "$r $g $b ")
                end
                println(io)
            end
        end
    end
end

end # module
