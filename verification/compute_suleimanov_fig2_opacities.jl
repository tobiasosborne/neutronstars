#!/usr/bin/env julia

# Compute the opacity components plotted in Suleimanov, Potekhin & Werner
# (2009) Fig. 2 for the paper's stated plasma parameters.

using NeutronStar
using NeutronStar.PhysicalConstants: h, keV, m_e, m_p
using NeutronStar.DielectricTensor: polarization_weights_full

const ROOT = dirname(dirname(@__FILE__))
const OUT = joinpath(ROOT, "verification", "data", "suleimanov_2009_fig2",
                     "computed_opacities.csv")

const T = 7.0e6
const RHO = 30.0
const B = 1.0e14
const M_H = m_p + m_e

function weighted_component(mode::Int, energy_keV::Float64, theta_deg::Float64,
                            component::Symbol)::Float64
    nu = energy_keV * keV / h
    omega = 2π * nu
    theta = deg2rad(theta_deg)
    n_e = RHO / M_H

    weights = polarization_weights_full(omega, B, theta, n_e)
    w = mode == 1 ? weights[1] : weights[2]

    sigma_sum = 0.0
    for (idx, alpha) in enumerate((-1, 0, 1))
        sigma = if component == :ff
            sigma_ff_alpha(alpha, omega, B, T, RHO)
        elseif component == :pp
            sigma_pp_alpha(alpha, omega, B, T, RHO)
        elseif component == :scattering
            sigma_scat_alpha(alpha, omega, B)
        else
            error("unknown component $component")
        end
        sigma_sum += w[idx] * sigma
    end

    return sigma_sum / M_H
end

function main()
    mkpath(dirname(OUT))
    energies = 10.0 .^ range(log10(0.01), log10(45.0), length=700)

    open(OUT, "w") do io
        println(io, "theta_deg,mode,energy_keV,ff_cm2_g,pp_cm2_g,absorption_cm2_g,scattering_cm2_g,total_cm2_g")
        for theta in (5.0, 45.0, 90.0)
            for mode in (1, 2)
                for energy in energies
                    ff = weighted_component(mode, energy, theta, :ff)
                    pp = weighted_component(mode, energy, theta, :pp)
                    scattering = weighted_component(mode, energy, theta, :scattering)
                    absorption = ff + pp
                    total = absorption + scattering
                    println(io, join((theta, mode, energy, ff, pp, absorption,
                                      scattering, total), ","))
                end
            end
        end
    end

    println("wrote ", relpath(OUT, ROOT))
end

main()
