#=
EOS and TOV verification against Potekhin et al. (2013).
Source: A&A 560, A48, Section 6 and Fig. 8.
=#

using NeutronStar
using Printf

function run_eos_verification()
    println("=" ^ 60)
    println("EOS & TOV VERIFICATION — Potekhin et al. (2013)")
    println("=" ^ 60)

    # === Test 1: EOS monotonicity and range ===
    println("\n--- Test 1: EOS monotonicity ---")
    P_prev = 0.0
    monotonic = true
    for logρ in 6.0:0.1:15.5
        P = pressure_of_density(10.0^logρ, BSk21_params)
        if P <= P_prev
            monotonic = false
            println("  FAIL: P not monotonic at log ρ = $logρ")
        end
        P_prev = P
    end
    println("  BSk21 P(ρ) monotonic: $(monotonic ? "PASS" : "FAIL")")

    # === Test 2: BSk21 M_max ===
    println("\n--- Test 2: BSk21 maximum mass ---")
    M_max = 0.0
    ρ_c_max = 0.0
    for logρc in 14.8:0.02:15.6
        ρc = 10.0^logρc
        try
            r = solve_tov(ρc, BSk21_params)
            if r.M_solar > M_max
                M_max = r.M_solar
                ρ_c_max = ρc
            end
        catch; end
    end
    M_max_expected = 2.27
    err = abs(M_max - M_max_expected) / M_max_expected * 100
    @printf("  M_max = %.4f M_☉ (expected: %.2f M_☉, error: %.2f%%)\n",
            M_max, M_max_expected, err)
    @printf("  ρ_c(M_max) = %.3e g/cm³ (expected: ~2.27e15)\n", ρ_c_max)
    println("  ", err < 2.0 ? "PASS" : "FAIL", " (tolerance: 2%)")

    # === Test 3: BSk20 M_max ===
    println("\n--- Test 3: BSk20 maximum mass ---")
    M_max_20 = 0.0
    for logρc in 14.8:0.02:15.6
        ρc = 10.0^logρc
        try
            r = solve_tov(ρc, BSk20_params)
            if r.M_solar > M_max_20
                M_max_20 = r.M_solar
            end
        catch; end
    end
    M_max_20_expected = 2.16
    err_20 = abs(M_max_20 - M_max_20_expected) / M_max_20_expected * 100
    @printf("  M_max = %.4f M_☉ (expected: %.2f M_☉, error: %.2f%%)\n",
            M_max_20, M_max_20_expected, err_20)
    println("  ", err_20 < 2.0 ? "PASS" : "FAIL")

    # === Test 4: BSk19 M_max ===
    println("\n--- Test 4: BSk19 maximum mass ---")
    M_max_19 = 0.0
    for logρc in 14.8:0.02:15.6
        ρc = 10.0^logρc
        try
            r = solve_tov(ρc, BSk19_params)
            if r.M_solar > M_max_19
                M_max_19 = r.M_solar
            end
        catch; end
    end
    M_max_19_expected = 1.86
    err_19 = abs(M_max_19 - M_max_19_expected) / M_max_19_expected * 100
    @printf("  M_max = %.4f M_☉ (expected: %.2f M_☉, error: %.2f%%)\n",
            M_max_19, M_max_19_expected, err_19)
    println("  ", err_19 < 2.0 ? "PASS" : "FAIL")

    # === Test 5: Canonical 1.4 M_☉ NS ===
    println("\n--- Test 5: Canonical 1.4 M_☉ NS (BSk21) ---")
    ρ_lo = 3e14
    ρ_hi = 1e15
    for _ in 1:50
        ρ_mid = sqrt(ρ_lo * ρ_hi)
        r = solve_tov(ρ_mid, BSk21_params)
        if r.M_solar < 1.4
            ρ_lo = ρ_mid
        else
            ρ_hi = ρ_mid
        end
    end
    ρ_14 = sqrt(ρ_lo * ρ_hi)
    r_14 = solve_tov(ρ_14, BSk21_params)
    @printf("  M = %.4f M_☉, R = %.2f km, ρ_c = %.3e g/cm³\n",
            r_14.M_solar, r_14.R_km, r_14.ρ_c)
    @printf("  g_surface = %.3e cm/s², compactness u = %.4f\n",
            r_14.g_surface, r_14.compactness)
    R_14_expected = 12.6
    err_R = abs(r_14.R_km - R_14_expected) / R_14_expected * 100
    @printf("  R error vs expected ~%.1f km: %.1f%%\n", R_14_expected, err_R)
    println("  ", err_R < 5.0 ? "PASS" : "FAIL", " (tolerance: 5%)")

    # === Test 6: Newtonian limit ===
    println("\n--- Test 6: Low-density / Newtonian limit ---")
    r_low = solve_tov(5e13, BSk21_params)
    @printf("  ρ_c = 5e13: M = %.4f M_☉, R = %.2f km, u = %.6f\n",
            r_low.M_solar, r_low.R_km, r_low.compactness)
    println("  Compactness u << 1: ", r_low.compactness < 0.01 ? "PASS" : "FAIL")

    println("\n" * "=" ^ 60)
    println("VERIFICATION COMPLETE")
    println("=" ^ 60)
end

run_eos_verification()
