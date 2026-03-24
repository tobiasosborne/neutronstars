#=
Physical constants in CGS units.
Source: CODATA 2018 (Tiesinga et al., Rev. Mod. Phys. 93, 025010, 2021).
All values are exact or given to full recommended precision.
=#

module PhysicalConstants

# Fundamental constants (CGS)
"Speed of light [cm/s]. CODATA 2018 exact."
const c = 2.99792458e10

"Gravitational constant [cm³ g⁻¹ s⁻²]. CODATA 2018."
const G = 6.67430e-8

"Planck constant [erg s]. CODATA 2018 exact."
const h = 6.62607015e-27

"Reduced Planck constant [erg s]. h/(2π)."
const ħ = h / (2π)

"Boltzmann constant [erg K⁻¹]. CODATA 2018 exact."
const k_B = 1.380649e-16

"Stefan-Boltzmann constant [erg cm⁻² s⁻¹ K⁻⁴]. σ = 2π⁵k⁴/(15h³c²)."
const σ_SB = 5.670374419e-5

"Electron mass [g]. CODATA 2018."
const m_e = 9.1093837015e-28

"Proton mass [g]. CODATA 2018."
const m_p = 1.67262192369e-24

"Atomic mass unit [g]. CODATA 2018."
const m_u = 1.66053906660e-24

"Elementary charge [esu = statcoulomb]. CODATA 2018."
const e_charge = 4.80320451e-10

"Thomson cross-section [cm²]. σ_T = 8πe⁴/(3m_e²c⁴). CODATA 2018."
const σ_T = 6.6524587321e-25

"Fine structure constant [dimensionless]. CODATA 2018."
const α_fine = 7.2973525693e-3

"Classical electron radius [cm]. r_e = e²/(m_e c²)."
const r_e = 2.8179403262e-13

# Astrophysical constants
"Solar mass [g]. IAU 2015 nominal."
const M_sun = 1.98892e33

"Parsec [cm]. IAU 2012 exact."
const pc = 3.0856775814913673e18

"Electron volt [erg]. CODATA 2018 exact."
const eV = 1.602176634e-12

"keV [erg]."
const keV = 1.0e3 * eV

# Derived constants for NS physics
"Schwarzschild radius per solar mass [cm]. r_s = 2GM/c²."
const r_s_per_Msun = 2G * M_sun / c^2  # ≈ 2.953 km

"Nuclear saturation density [g cm⁻³]. ρ_s = n_s × m_u."
const ρ_nuc = 2.8e14

"Critical magnetic field (QED) [G]. B_Q = m_e²c³/(eℏ)."
const B_Q = m_e^2 * c^3 / (e_charge * ħ)  # ≈ 4.414e13 G

"Cyclotron energy per field [eV/G]. ℏω_c/B = eℏ/(m_e c)."
const ħωc_per_B = e_charge * ħ / (m_e * c) / eV  # eV per Gauss

end # module
