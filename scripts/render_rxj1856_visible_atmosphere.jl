using Printf
using NeutronStar
using NeutronStar.PhysicalConstants: G, c, M_sun
using NeutronStar.GauntFactor: load_gaunt_table

const OUT_PREFIX = "output/rxj1856_visible_atmosphere"

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

function render_rxj1856_visible(; image_n::Int=256,
                                atmosphere_k::Int=80,
                                atmosphere_m::Int=8,
                                atmosphere_n::Int=200,
                                visible_nu_n::Int=180)
    mkpath(dirname(OUT_PREFIX))

    # RX J1856-like simple model from Ho et al. 2007:
    # T_inf = 4.3e5 K, R_inf = 17 km, z = 0.22.
    # The current B>0 atmosphere is not flux-normalized, so this render uses
    # the working non-magnetic hydrogen atmosphere at the corresponding surface
    # T_eff = T_inf * (1 + z).
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

    @printf("RX J1856 visible atmosphere render\n")
    @printf("  M=%.3f Msun, R=%.3f km, z=%.3f, T_eff=%.3e K, g_s=%.3e\n",
            M_solar, R_km, z_gr, T_eff, g_s)

    gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
    atmosphere = solve_atmosphere(T_eff, g_s, gaunt;
                                  nν=atmosphere_k,
                                  M=atmosphere_m,
                                  N=atmosphere_n,
                                  max_iter=30,
                                  tol=1e-4,
                                  verbose=true)

    image = trace_image(M, R, 0.0, image_n)

    lambda_nm = range(830.0, 360.0, length=visible_nu_n)
    nu_grid = [c / (lambda * 1e-7) for lambda in lambda_nm]

    spectra = Array{Vector{Float64}}(undef, image_n, image_n)
    Y_values = Float64[]
    cmfs = load_cie_cmfs(joinpath(dirname(pathof(NeutronStar)), "data", "cvrl_cie1931_2deg.csv"))

    for j in 1:image_n, i in 1:image_n
        ray = image.rays[i, j]
        if !ray.hits
            spectra[i, j] = zeros(visible_nu_n)
            continue
        end

        mu = clamp(ray.cos_α, 1e-6, 1.0)
        surface_spectrum = rt_emergent_spectrum(atmosphere, mu, nu_grid .* redshift)
        spectra[i, j] = surface_spectrum ./ redshift^3

        _, Y, _ = spectrum_to_XYZ(nu_grid, spectra[i, j], cmfs)
        push!(Y_values, Y)
    end

    L_avg = isempty(Y_values) ? 1.0 : exp(sum(log.(max.(Y_values, 1e-300))) / length(Y_values))
    rgb = Array{NTuple{3,Int}}(undef, image_n, image_n)

    for j in 1:image_n, i in 1:image_n
        if !image.rays[i, j].hits
            rgb[i, j] = (0, 0, 0)
        else
            rgb[i, j] = spectrum_to_sRGB(nu_grid, spectra[i, j], cmfs;
                                         exposure=1.0,
                                         L_avg=L_avg,
                                         key=0.36)
        end
    end

    ppm_path = OUT_PREFIX * ".ppm"
    save_ppm(ppm_path, rgb)

    open(OUT_PREFIX * ".txt", "w") do io
        println(io, "RX J1856 visible atmosphere render")
        println(io, "object = RX J1856.5-3754 / RX J185635-3754")
        println(io, "atmosphere = non-magnetic fully ionized H, solve_atmosphere")
        println(io, "reason_B0 = repo notes/failures.md records B>0 atmosphere flux not normalized")
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
    return ppm_path
end

render_rxj1856_visible()
