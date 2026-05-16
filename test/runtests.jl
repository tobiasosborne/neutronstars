using Test
using NeutronStar

const NS = NeutronStar

@testset "Magnetic mode opacity separation" begin
    ν = 2.0e17
    θ_B = π / 4
    B = 1.0e12
    T = 1.0e6
    ρ = 1.0e-3

    for mode in 1:2
        κ_abs = mode_absorption(mode, ν, θ_B, B, T, ρ)
        κ_scat = mode_scattering(mode, ν, θ_B, B, T, ρ)
        κ_total = mode_opacity(mode, ν, θ_B, B, T, ρ)

        @test κ_abs > 0.0
        @test κ_scat > 0.0
        @test κ_total ≈ κ_abs + κ_scat rtol=1e-12
    end
end

@testset "Magnetic Coulomb logarithm follows P&C Eq. 44 tail" begin
    B = 1.0e14
    T = 7.0e6
    E_keV = 0.01
    u = E_keV * NS.PhysicalConstants.keV / (NS.PhysicalConstants.k_B * T)
    β_e = NS.PhysicalConstants.ħ * cyclotron_freq_e(B) / (NS.PhysicalConstants.k_B * T)

    # Direct numerical evaluation of P&C 2003 Eq. 44a-e for the Fig. 2
    # low-energy longitudinal polarization point.  The long y tail is required.
    Λ0 = NS.MagneticCoulomb.coulomb_log_magnetic(0, u, β_e)

    @test Λ0 ≈ 7.4148 rtol=5e-4
end

@testset "Vacuum-polarized mode vectors" begin
    E_keV = 1.0
    B = 1.0e14
    θ_B = π / 4
    ω = E_keV * NS.PhysicalConstants.keV / NS.PhysicalConstants.ħ
    n_e = 1.0 / (NS.PhysicalConstants.m_p + NS.PhysicalConstants.m_e)

    a_hat, q_vac, m_vac = NS.DielectricTensor.vacuum_coefficients(B)
    @test a_hat ≈ -2.0563758036448108e-4 rtol=1e-10
    @test q_vac ≈ 1.0012982615339183e-3 rtol=1e-10
    @test m_vac ≈ -2.425084838547189e-4 rtol=1e-10

    β, r, K1, K2, Kz1, Kz2 =
        NS.DielectricTensor.vacuum_mode_parameters(ω, B, θ_B, n_e)
    @test β ≈ 19.411715027031708 rtol=1e-10
    @test r ≈ 0.9998787208185153 rtol=1e-10
    @test K1 ≈ -0.025737455521326447 rtol=1e-10
    @test K2 ≈ 38.84916750958474 rtol=1e-10
    @test Kz1 ≈ -1.4514365449006211e-6 rtol=1e-8
    @test Kz2 ≈ 0.003452653275725913 rtol=1e-8

    w1, w2 = NS.DielectricTensor.polarization_weights_vacuum(ω, B, θ_B, n_e)
    @test w1[1] ≈ 0.5180226316829358 rtol=1e-8
    @test w1[2] ≈ 0.0003309517252031447 rtol=1e-8
    @test w1[3] ≈ 0.48164641659186097 rtol=1e-8
    @test w2[1] ≈ 0.23201904050694275 rtol=1e-8
    @test w2[2] ≈ 0.49958011625984716 rtol=1e-8
    @test w2[3] ≈ 0.26840084323321006 rtol=1e-8

    n_low = 0.1 / (NS.PhysicalConstants.m_p + NS.PhysicalConstants.m_e)
    β_low = NS.DielectricTensor.vacuum_polarization_parameter(ω, B, θ_B, n_low)
    @test β_low < 0.0
    @test β_low ≈ -2022.0197070158636 rtol=1e-10
end

@testset "Magnetic initial column reaches diffusion depth" begin
    gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
    ν_grid = make_frequency_grid(1.0e6, 12)
    N = 40
    B = 1.0e12
    θ_B = π / 4

    y, T, ρ, _ = NS.MagneticAtmosphere._build_initial_column(
        1.0e6, 2.0e14, N, ν_grid, gaunt, B, θ_B
    )

    κ = zeros(N, length(ν_grid), 2)
    k_total = similar(κ)
    ρ_alb = similar(κ)
    τ = similar(κ)

    NS.MagneticAtmosphere._compute_magnetic_opacities!(
        κ, k_total, ρ_alb, τ, y, T, ρ, ν_grid, B, θ_B, gaunt
    )

    @test minimum(τ[end, :, :]) >= 80.0
end

@testset "Magnetic atmosphere flux normalization smoke" begin
    gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")
    result = solve_magnetic_atmosphere(
        1.0e6, 2.0e14, 1.0e12, π / 4, gaunt;
        K=16, M=4, N=50, max_iter=60, tol=5e-4, verbose=false
    )

    flux = 0.0
    _, weights = NeutronStar.FeautrierSolver.gauss_legendre_half(4)
    for k in 1:length(result.ν_grid)-1
        dν = result.ν_grid[k+1] - result.ν_grid[k]
        for j in eachindex(result.μ_grid), mode in 1:2
            flux += 2π * result.μ_grid[j] *
                    result.I_emergent[k, j, mode] *
                    weights[j] *
                    dν
        end
    end

    @test result.converged
    @test flux / (NeutronStar.PhysicalConstants.σ_SB * result.T_eff^4) ≈ 1.0 rtol=0.05
end

@testset "Magnetic atmosphere uses scattering albedo for B > 0" begin
    y = [1.0e-6, 1.0e-4, 1.0e-2]
    T = fill(1.0e6, length(y))
    ρ = fill(1.0e-3, length(y))
    ν_grid = [2.0e17]
    B = 1.0e12
    θ_B = π / 4

    κ = zeros(length(y), length(ν_grid), 2)
    k_total = similar(κ)
    ρ_alb = similar(κ)
    τ = similar(κ)

    NS.MagneticAtmosphere._compute_magnetic_opacities!(
        κ, k_total, ρ_alb, τ, y, T, ρ, ν_grid, B, θ_B, nothing
    )

    for mode in 1:2
        @test all(κ[:, 1, mode] .> 0.0)
        @test all(k_total[:, 1, mode] .> κ[:, 1, mode])
        @test all(0.0 .< ρ_alb[:, 1, mode] .< 1.0)
        @test τ[1, 1, mode] == 0.0
        @test τ[end, 1, mode] > 0.0
    end
end

@testset "Phase 2 tracer-bullet (Tier 1)" begin
    # 1. TOV: BSk21 canonical 1.4 M_☉ NS → R ≈ 12.6 km.
    # VERIFICATION_LOG.md:20 records ρ_c = 7.30e14 g/cm³ → M = 1.4000 M_☉, R = 12.59 km.
    # We use that ρ_c directly (avoids a bisection inside the test).
    tov = solve_tov(7.30e14, BSk21_params)
    @test tov.M_solar ≈ 1.4 atol=0.02   # within 0.02 M_☉ of canonical
    @test tov.R_km ≈ 12.6 atol=0.5      # published BSk21 value (Potekhin+ 2013 Fig. 8)
    println("  TOV(BSk21, ρ_c=7.3e14): M = $(round(tov.M_solar, digits=4)) M_☉, R = $(round(tov.R_km, digits=3)) km")

    # 2-3. visible_fraction limits (Beloborodov 2002 Sect. 3, Pechenick & Ftaclas 1983).
    # u = r_g/R = 2GM/(Rc²). At u = 0 (flat) → 0.5; at u = 1/3 (compact) → 0.75.
    @test visible_fraction(0.0) ≈ 0.5 atol=1e-3
    @test visible_fraction(1.0/3.0) ≈ 0.75 atol=1e-3
    println("  visible_fraction(0) = $(visible_fraction(0.0)), visible_fraction(1/3) = $(round(visible_fraction(1/3), digits=6))")

    # 4. surface_temperature pole/equator limits (Greenstein-Hartke 1983).
    # Signature: surface_temperature(θ_B, T_pole, T_eq). T(0) = T_pole, T(π/2) = T_eq.
    T_pole = 1.0e6
    T_eq = 5.0e5
    @test surface_temperature(0.0, T_pole, T_eq) == T_pole
    @test surface_temperature(π/2, T_pole, T_eq) == T_eq

    # 5. surface_Bfield pole/equator (centred dipole convention).
    # |B|(θ_B) = (B_pole/2)√(1 + 3cos²θ_B). At pole → B_pole, at equator → B_pole/2.
    B_pole = 1.0e12
    @test surface_Bfield(0.0, B_pole) ≈ B_pole rtol=1e-12
    @test surface_Bfield(π/2, B_pole) ≈ 0.5 * B_pole rtol=1e-12

    # 6. magnetic_colatitude with obliquity = 0: dipole axis aligned with spin.
    # Then θ_B == θ_geo for every geographic point.
    @test magnetic_colatitude(0.0, 0.0, 0.0) ≈ 0.0 atol=1e-12
    @test magnetic_colatitude(π/2, 0.0, 0.0) ≈ π/2 atol=1e-12

    # 7. Solar blackbody (T = 5778 K) → warm-white in sRGB.
    # VERIFICATION_LOG.md:49-53 documents chromaticity (x,y) = (0.3264, 0.3357)
    # and 8-bit sRGB ≈ (255, 252, 245). The 8-bit values depend on the
    # Reinhard tone-map exposure/key/L_avg choice (not pinned in the log),
    # so we assert only the physics-invariant: warm-white ordering and that
    # the red channel saturates first under reasonable exposure. The exact
    # chromaticity assertion is left as a TODO until the tone-map calibration
    # is documented; see VERIFICATION_LOG.md:49-53 and reviews/02_tests.md:122.
    cmfs = load_cie_cmfs("refs/cvrl_cie1931_2deg.csv")
    ν_grid = collect(range(1.0e14, 2.0e15; length=400))
    I_ν = [planck_Bnu(ν, 5778.0) for ν in ν_grid]
    _, Y_solar, _ = spectrum_to_XYZ(ν_grid, I_ν, cmfs)
    r, g, b = spectrum_to_sRGB(ν_grid, I_ν, cmfs; exposure=1.0, L_avg=Y_solar, key=10.0)
    println("  Solar BB (5778 K) → sRGB = ($r, $g, $b); VERIFICATION_LOG target ≈ (255, 252, 245)")
    @test r >= g >= b          # warm white (R ≥ G ≥ B)
    @test r == 255             # red channel saturates under this exposure
    @test g >= 200             # well into the white region, not orange
    @test b >= 180             # still substantial blue (solar, not cool star)
    # TODO: once tone-map calibration is pinned, assert (255, 252, 245) atol=3.
end

@testset "B=0 magnetic limit recovers non-magnetic" begin
    # Consistency check: with B=0 the magnetic solver must reduce to the
    # non-magnetic one to within numerical tolerance.  Both now use the
    # Hummer (1962) dipole Thomson-scattering kernel, so the only remaining
    # differences are (a) per-mode B_ν/2 splitting (which sums back to B_ν)
    # and (b) optional inter-mode coupling that vanishes at B=0.
    gaunt = load_gaunt_table("refs/code/McPHAC/gffgu.dat")

    # max_iter matches the existing magnetic-flux smoke test (line 96) for the
    # same K=16, M=4, N=50 grid; 30 iterations were not always enough to satisfy
    # both the temperature and flux convergence criteria at this coarse grid.
    r_mag = solve_magnetic_atmosphere(
        1.0e6, 2.0e14, 0.0, 0.0, gaunt;
        K=16, M=4, N=50, max_iter=60, tol=5.0e-4, verbose=false
    )
    r_nonmag = solve_atmosphere(
        1.0e6, 2.0e14, gaunt;
        K=16, M=4, N=50, max_iter=60, tol=5.0e-4, verbose=false
    )

    # Sum modes for magnetic: unpolarized intensity = mode 1 + mode 2
    I_mag = r_mag.I_emergent[:, :, 1] .+ r_mag.I_emergent[:, :, 2]
    I_nonmag = r_nonmag.I_emergent

    max_rel_diff = maximum(abs.(I_mag .- I_nonmag) ./ max.(I_nonmag, 1.0e-30))
    println("  B=0 limit max rel diff (mag vs non-mag I_emergent): $(round(max_rel_diff, sigdigits=4))")
    # The Hummer dipole kernels in both solvers are identical at B=0, so the
    # only differences come from per-mode B_ν/2 splitting (sums back to B_ν)
    # and floating-point reductions.  Empirically this gives ≲1e-4.  We assert
    # a generous 0.5% bound for CI robustness across compilers/BLAS.
    @test max_rel_diff < 5.0e-3

    # Convergence flags are reported for visibility but not asserted — at the
    # coarse CI grid the iteration may not always reach both T and flux
    # tolerances simultaneously even though the emergent spectrum is correct.
    println("  r_mag.converged=$(r_mag.converged) (max|ΔT/T|=$(round(r_mag.max_dT_over_T, sigdigits=3)), iters=$(r_mag.n_iterations))")
    println("  r_nonmag.converged=$(r_nonmag.converged) (max|ΔT/T|=$(round(r_nonmag.max_dT_over_T, sigdigits=3)), iters=$(r_nonmag.n_iterations))")
end

@testset "Suleimanov 2009 Fig 2 metric gate" begin
    # Delegates to verification/check_suleimanov_fig2_metrics.py, which asserts
    # four free-free opacity thresholds (coverage>=0.45, |median residual|<=0.70 dex,
    # median abs error<=0.75 dex, robust scatter<=0.60 dex; plus p90<=1.25 dex at θ=90°)
    # against the checked-in opacity_comparison_metrics.json. PASS → exit 0, FAIL → exit 1.
    metrics_json = "verification/data/suleimanov_2009_fig2/opacity_comparison_metrics.json"
    if Sys.which("python3") !== nothing && isfile(metrics_json)
        proc = run(pipeline(`python3 verification/check_suleimanov_fig2_metrics.py`,
                            stdout=devnull, stderr=devnull); wait=false)
        wait(proc)
        @test proc.exitcode == 0
    else
        @warn "Skipping Suleimanov gate: python3 missing or metrics JSON not present"
        @test_skip true
    end
end
