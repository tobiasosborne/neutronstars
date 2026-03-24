#=
Comparison harness: our computed Rosseland opacities vs Potekhin tables.

Currently only compares columns 13-14 (lg K_∥, lg K_⊥) since the EOS
modules (cols 2-12) haven't been implemented yet.

Usage:
  julia --project=. verification/potekhin_table_comparison.jl
=#

include(joinpath(@__DIR__, "..", "src", "eos", "potekhin_table_reader.jl"))
using .PotekhinTableReader

using NeutronStar
using NeutronStar.PhysicalConstants: k_B, h, m_e, m_p
using Printf

const m_H = m_p + m_e

"""
    compare_opacity(lgB, lgT, lgR; N_ν=60) → (lgK0_ours, lgK1_ours)

Compute our Rosseland mean opacities for one (B, T, ρ) point.
"""
function compare_opacity(lgB::Float64, lgT::Float64, lgR::Float64; N_ν::Int=60)
    B = 10.0^lgB
    T = 10.0^lgT
    T6 = T / 1e6
    ρ = 10.0^lgR * T6^3

    # Construct frequency grid
    ν_min = 0.05 * k_B * T / h
    ν_max = 120.0 * k_B * T / h
    ν_grid = [10.0^logν for logν in range(log10(ν_min), log10(ν_max), length=N_ν)]

    K_par, K_perp = rosseland_magnetic(B, T, ρ, ν_grid)

    lgK0 = K_par > 0 ? log10(K_par) : -30.0
    lgK1 = K_perp > 0 ? log10(K_perp) : -30.0

    return lgK0, lgK1
end

"""
    run_comparison(; test_lgBs, N_ν, max_points_per_iso)

Compare our opacities against Potekhin tables for selected B values.
"""
function run_comparison(;
    test_lgBs::Vector{Float64} = [10.5, 11.0, 12.0, 13.0, 13.5, 14.0, 15.0],
    N_ν::Int = 60,
    max_points_per_iso::Int = 10  # subsample for speed
)
    println("Loading Potekhin tables...")
    all_tables = load_all_tables(joinpath(@__DIR__, "..", "refs", "potekhin_tables"))
    merged = merge_tables(all_tables)
    println("Loaded $(length(merged)) lgB values\n")

    println("="^100)
    @printf("%-6s %-6s %-6s │ %8s %8s %6s │ %8s %8s %6s │ %6s\n",
            "lgB", "lgT", "lgR", "K∥_ours", "K∥_tab", "Δ∥", "K⊥_ours", "K⊥_tab", "Δ⊥", "x_H")
    println("-"^100)

    # Collect error statistics
    errors_par = Float64[]
    errors_perp = Float64[]
    errors_par_ionized = Float64[]
    errors_perp_ionized = Float64[]

    for lgB_target in test_lgBs
        # Find closest lgB in merged tables
        lgB_key = nothing
        for k in keys(merged)
            if abs(k - lgB_target) < 0.05
                lgB_key = k
                break
            end
        end
        lgB_key === nothing && continue

        table = merged[lgB_key]

        # Pick a few representative isotherms
        n_iso = length(table.isotherms)
        iso_indices = unique(clamp.(round.(Int, range(1, n_iso, length=min(3, n_iso))), 1, n_iso))

        for idx in iso_indices
            iso = table.isotherms[idx]
            lgT = iso.lgT
            lgB = iso.lgB

            # Subsample density points
            n_entries = length(iso.entries)
            step = max(1, n_entries ÷ max_points_per_iso)
            sample_indices = 1:step:n_entries

            for i in sample_indices
                e = iso.entries[i]
                lgR = e.lgR

                lgK0_ours, lgK1_ours = try
                    compare_opacity(lgB, lgT, lgR; N_ν=N_ν)
                catch ex
                    @warn "Failed at lgB=$lgB lgT=$lgT lgR=$lgR: $ex"
                    (-30.0, -30.0)
                end

                # Skip comparison where table has effectively zero opacity
                if e.lgK0 < -9.0 && e.lgK1 < -9.0
                    continue
                end

                Δ_par = lgK0_ours - e.lgK0
                Δ_perp = lgK1_ours - e.lgK1

                is_ionized = e.x_H < 0.01

                @printf("%-6.1f %-6.3f %+5.1f  │ %8.3f %8.3f %+6.2f │ %8.3f %8.3f %+6.2f │ %6.3f %s\n",
                        lgB, lgT, lgR, lgK0_ours, e.lgK0, Δ_par,
                        lgK1_ours, e.lgK1, Δ_perp, e.x_H,
                        is_ionized ? "" : " [partial]")

                push!(errors_par, abs(Δ_par))
                push!(errors_perp, abs(Δ_perp))
                if is_ionized
                    push!(errors_par_ionized, abs(Δ_par))
                    push!(errors_perp_ionized, abs(Δ_perp))
                end
            end
        end
    end

    println("="^100)
    println("\nError summary (absolute error in dex):")
    println("-"^60)

    function stats(v, label)
        isempty(v) && return
        @printf("  %-30s  N=%4d  med=%.2f  mean=%.2f  max=%.2f\n",
                label, length(v), median_sorted(sort(v)),
                sum(v)/length(v), maximum(v))
    end

    stats(errors_par, "K∥ (all points)")
    stats(errors_perp, "K⊥ (all points)")
    stats(errors_par_ionized, "K∥ (fully ionized, x_H<0.01)")
    stats(errors_perp_ionized, "K⊥ (fully ionized, x_H<0.01)")
end

function median_sorted(v)
    n = length(v)
    n == 0 && return NaN
    n % 2 == 1 ? v[(n+1)÷2] : (v[n÷2] + v[n÷2+1]) / 2
end

# Run
run_comparison()
