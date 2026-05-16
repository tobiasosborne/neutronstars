using Printf
using NeutronStar
using NeutronStar.PhysicalConstants: G, c, M_sun
using NeutronStar.GauntFactor: load_gaunt_table

const OUT_PREFIX = "output/rxj1856_visible_magnetic_atmosphere"

function save_ppm(path::String, rgb)
    n_x, n_y = size(rgb)
    open(path, "w") do io
        println(io, "P3")
        println(io, "$n_x $n_y")
        println(io, "255")
        for j in n_y:-1:1
            for i in 1:n_x
                r, g, b = rgb[i, j]
                print(io, "$r $g $b ")
            end
            println(io)
        end
    end
end

function magnetic_spectrum(result, cos_θe, ν_out)
    μ = result.μ_grid
    M = length(μ)
    j_lo = 1
    for j in 1:M-1
        if μ[j] <= cos_θe <= μ[j+1]
            j_lo = j
            break
        end
    end
    if cos_θe < μ[1]; j_lo = 1; end
    if cos_θe > μ[end]; j_lo = M - 1; end
    j_hi = min(j_lo + 1, M)
    t_μ = clamp((cos_θe - μ[j_lo]) / (μ[j_hi] - μ[j_lo]), 0.0, 1.0)

    ν_in = result.ν_grid
    K = length(ν_in)
    I_out = zeros(length(ν_out))

    for i in eachindex(ν_out)
        ν = ν_out[i]
        if ν <= ν_in[1]
            I_lo = sum(result.I_emergent[1, j_lo, mode] for mode in 1:2)
            I_hi = sum(result.I_emergent[1, j_hi, mode] for mode in 1:2)
            I_out[i] = ((1.0 - t_μ) * I_lo + t_μ * I_hi) * (ν / ν_in[1])^2
        elseif ν >= ν_in[end]
            I_out[i] = 0.0
        else
            k_lo = clamp(searchsortedlast(ν_in, ν), 1, K-1)
            k_hi = k_lo + 1
            t_ν = (log(ν) - log(ν_in[k_lo])) / (log(ν_in[k_hi]) - log(ν_in[k_lo]))

            I_ll = sum(result.I_emergent[k_lo, j_lo, mode] for mode in 1:2)
            I_hl = sum(result.I_emergent[k_hi, j_lo, mode] for mode in 1:2)
            I_lh = sum(result.I_emergent[k_lo, j_hi, mode] for mode in 1:2)
            I_hh = sum(result.I_emergent[k_hi, j_hi, mode] for mode in 1:2)

            I_lo = I_ll > 0 && I_hl > 0 ? exp((1.0-t_ν)*log(I_ll) + t_ν*log(I_hl)) :
                   (1.0-t_ν)*I_ll + t_ν*I_hl
            I_hi = I_lh > 0 && I_hh > 0 ? exp((1.0-t_ν)*log(I_lh) + t_ν*log(I_hh)) :
                   (1.0-t_ν)*I_lh + t_ν*I_hh
            I_out[i] = (1.0 - t_μ) * I_lo + t_μ * I_hi
        end
    end

    return I_out
end

function render_rxj1856_visible_magnetic(; image_n::Int=256)
    mkpath(dirname(OUT_PREFIX))

    z_gr = 0.22
    redshift = 1.0 + z_gr
    T_inf = 4.3e5
    T_eff = T_inf * redshift
    R_inf_km = 17.0
    R_km = R_inf_km / redshift
    compactness = 1.0 - redshift^(-2)
    r_g_per_msun_km = 2.0 * G * M_sun / c^2 / 1e5
    M_solar = compactness * R_km / r_g_per_msun_km
    M = M_solar * M_sun
    R = R_km * 1e5
    g_s = G * M / R^2 * redshift

    # This is a verified-convergent magnetic case.  The stronger RX J1856
    # inferred field (~3e12 G) is flux-normalized by the current code but does
    # not yet satisfy the local temperature-correction tolerance.
    B = 1.0e12
    θ_B = π / 4

    @printf("RX J1856 visible magnetic atmosphere render\n")
    @printf("  M=%.3f Msun, R=%.3f km, z=%.3f, T_eff=%.3e K, g_s=%.3e, B=%.2e G\n",
            M_solar, R_km, z_gr, T_eff, g_s, B)

    gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
    atmosphere = solve_magnetic_atmosphere(T_eff, g_s, B, θ_B, gaunt;
                                           K=50, M=8, N=100,
                                           max_iter=45,
                                           tol=7e-4,
                                           flux_tol=1e-2,
                                           flux_damping=0.5,
                                           verbose=true)

    image = trace_image(M, R, 0.0, image_n)
    lambda_nm = range(830.0, 360.0, length=180)
    ν_grid = [c / (lambda * 1e-7) for lambda in lambda_nm]
    cmfs = load_cie_cmfs("refs/cvrl_cie1931_2deg.csv")

    spectra = Array{Vector{Float64}}(undef, image_n, image_n)
    Y_values = Float64[]
    for j in 1:image_n, i in 1:image_n
        ray = image.rays[i, j]
        if !ray.hits
            spectra[i, j] = zeros(length(ν_grid))
            continue
        end
        surface = magnetic_spectrum(atmosphere, clamp(ray.cos_α, 1e-6, 1.0), ν_grid .* redshift)
        spectra[i, j] = surface ./ redshift^3
        _, Y, _ = spectrum_to_XYZ(ν_grid, spectra[i, j], cmfs)
        push!(Y_values, Y)
    end

    L_avg = exp(sum(log.(max.(Y_values, 1e-300))) / length(Y_values))
    rgb = Array{NTuple{3,Int}}(undef, image_n, image_n)
    for j in 1:image_n, i in 1:image_n
        if !image.rays[i, j].hits
            rgb[i, j] = (0, 0, 0)
        else
            rgb[i, j] = spectrum_to_sRGB(ν_grid, spectra[i, j], cmfs;
                                         exposure=1.0, L_avg=L_avg, key=0.36)
        end
    end

    ppm_path = OUT_PREFIX * ".ppm"
    save_ppm(ppm_path, rgb)

    open(OUT_PREFIX * ".txt", "w") do io
        println(io, "RX J1856 visible magnetic atmosphere render")
        println(io, "object = RX J1856.5-3754 / RX J185635-3754")
        println(io, "atmosphere = magnetic fully ionized H, solve_magnetic_atmosphere")
        println(io, "B_G = $B")
        println(io, "theta_B_rad = $θ_B")
        println(io, "T_inf_K = $T_inf")
        println(io, "T_eff_surface_K = $T_eff")
        println(io, "z_gr = $z_gr")
        println(io, "R_inf_km = $R_inf_km")
        println(io, "R_km = $R_km")
        println(io, "M_solar = $M_solar")
        println(io, "g_s_cgs = $g_s")
        println(io, "image_n = $image_n")
        println(io, "visible_lambda_nm = 360..830")
        println(io, "atmosphere_converged = $(atmosphere.converged)")
        println(io, "atmosphere_iterations = $(atmosphere.n_iterations)")
        println(io, "atmosphere_max_dT_over_T = $(atmosphere.max_dT_over_T)")
        println(io, "visible_fraction = $(image.visible_frac)")
    end

    @printf("  wrote %s and %s\n", ppm_path, OUT_PREFIX * ".txt")
end

render_rxj1856_visible_magnetic()
