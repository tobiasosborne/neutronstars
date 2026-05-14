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
    @test flux / (NeutronStar.PhysicalConstants.σ_SB * result.T_eff^4) ≈ 1.0 rtol=0.25
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
