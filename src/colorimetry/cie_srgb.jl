#=
Colorimetry: CIE 1931 → XYZ → sRGB → tone mapping.

Sources:
- CIE 1931 2° observer: CVRL data (cvrl.ioo.ucl.ac.uk). Local: refs/cvrl_cie1931_2deg.csv
- XYZ → sRGB matrix: IEC 61966-2-1:1999 (exact values)
- sRGB gamma: IEC 61966-2-1:1999
- Tone mapping: Reinhard et al. SIGGRAPH 2002. Local: refs/reinhard_2002_tonemapping.pdf
=#

module CIE_sRGB

using ..PhysicalConstants: c, h

export spectrum_to_XYZ, XYZ_to_linear_sRGB, linear_sRGB_to_sRGB
export tone_map_reinhard, spectrum_to_sRGB
export load_cie_cmfs

"""
CIE 1931 2° colour matching functions, interpolated to a regular grid.
"""
struct CIE_CMFs
    λ::Vector{Float64}  # wavelengths [nm]
    x̄::Vector{Float64}
    ȳ::Vector{Float64}
    z̄::Vector{Float64}
end

"""
    load_cie_cmfs(path) → CIE_CMFs

Load CIE 1931 2° CMFs from CVRL CSV file.
Format: λ[nm], x̄, ȳ, z̄ (one row per nm, 360-830 nm).
"""
function load_cie_cmfs(path::String)::CIE_CMFs
    λ = Float64[]
    x̄ = Float64[]
    ȳ = Float64[]
    z̄ = Float64[]

    for line in eachline(path)
        parts = split(strip(line), ',')
        length(parts) >= 4 || continue
        push!(λ, parse(Float64, parts[1]))
        push!(x̄, parse(Float64, parts[2]))
        push!(ȳ, parse(Float64, parts[3]))
        push!(z̄, parse(Float64, parts[4]))
    end

    @assert length(λ) > 100 "Too few CMF data points: $(length(λ))"
    return CIE_CMFs(λ, x̄, ȳ, z̄)
end

"""
    spectrum_to_XYZ(ν_grid, I_ν, cmfs) → (X, Y, Z)

Integrate observed spectrum against CIE colour matching functions.

X = ∫ I_ν(ν) x̄(ν) dν  (and similarly Y, Z)

Converts from frequency domain (I_ν) to wavelength domain internally.
"""
function spectrum_to_XYZ(ν_grid::AbstractVector{Float64},
                          I_ν::AbstractVector{Float64},
                          cmfs::CIE_CMFs)
    @assert length(ν_grid) == length(I_ν)

    X = 0.0
    Y = 0.0
    Z = 0.0

    # Convert I_ν(ν) to I_λ(λ) using I_λ = I_ν × ν²/c = I_ν × c/λ²
    # Integrate in wavelength space where CMFs are defined

    for i in 1:length(cmfs.λ)-1
        λ_nm = cmfs.λ[i]
        λ_cm = λ_nm * 1e-7  # nm → cm
        ν = c / λ_cm  # Hz

        # Interpolate I_ν at this frequency
        I = _interp_spectrum(ν, ν_grid, I_ν)

        # I_λ = I_ν × c/λ² (change of variables dν = -c/λ² dλ)
        I_λ = I * c / λ_cm^2

        dλ = (cmfs.λ[i+1] - cmfs.λ[i]) * 1e-7  # nm → cm

        X += I_λ * cmfs.x̄[i] * dλ
        Y += I_λ * cmfs.ȳ[i] * dλ
        Z += I_λ * cmfs.z̄[i] * dλ
    end

    return X, Y, Z
end

"""
Linear interpolation of spectrum in log-frequency space.
"""
function _interp_spectrum(ν::Float64, ν_grid::AbstractVector{Float64},
                           I_ν::AbstractVector{Float64})::Float64
    # ν_grid should be monotonically increasing
    if ν <= ν_grid[1]
        return I_ν[1]
    elseif ν >= ν_grid[end]
        return I_ν[end]
    end

    # Binary search
    lo, hi = 1, length(ν_grid)
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        if ν_grid[mid] <= ν
            lo = mid
        else
            hi = mid
        end
    end

    # Linear interpolation in log space
    t = (log(ν) - log(ν_grid[lo])) / (log(ν_grid[hi]) - log(ν_grid[lo]))
    # Interpolate in log I to handle large dynamic range
    I_lo = max(I_ν[lo], 1e-300)
    I_hi = max(I_ν[hi], 1e-300)
    return exp(log(I_lo) + t * (log(I_hi) - log(I_lo)))
end

"""
    XYZ_to_linear_sRGB(X, Y, Z) → (R, G, B)

Convert CIE XYZ to linear sRGB using the IEC 61966-2-1 matrix.
Source: IEC 61966-2-1:1999.
"""
function XYZ_to_linear_sRGB(X::Float64, Y::Float64, Z::Float64)
    # Exact IEC 61966-2-1 matrix (D65 illuminant)
    R =  3.2406255 * X - 1.537208  * Y - 0.4986286 * Z
    G = -0.9689307 * X + 1.8757561 * Y + 0.0415175 * Z
    B =  0.0557101 * X - 0.2040211 * Y + 1.0569959 * Z
    return R, G, B
end

"""
    linear_sRGB_to_sRGB(C_lin) → C_srgb

Apply sRGB gamma curve. IEC 61966-2-1:1999.
"""
function linear_sRGB_to_sRGB(C_lin::Float64)::Float64
    if C_lin <= 0.0031308
        return 12.92 * C_lin
    else
        return 1.055 * C_lin^(1.0/2.4) - 0.055
    end
end

"""
    tone_map_reinhard(L; key=0.18) → L_display

Reinhard global tone mapping operator.
L_display = L / (1 + L)

After scaling: L_scaled = key × L / L_avg, then map.
Source: Reinhard et al. (2002) Eq. 1-2.
"""
function tone_map_reinhard(L::Float64; key::Float64=0.18,
                            L_avg::Float64=1.0)::Float64
    L_scaled = key * L / max(L_avg, 1e-300)
    return L_scaled / (1.0 + L_scaled)
end

"""
    spectrum_to_sRGB(ν_grid, I_ν, cmfs; exposure, tone_map) → (r, g, b) ∈ [0,255]

Full pipeline: spectrum → XYZ → linear sRGB → tone map → gamma → 8-bit.
"""
function spectrum_to_sRGB(ν_grid::AbstractVector{Float64},
                           I_ν::AbstractVector{Float64},
                           cmfs::CIE_CMFs;
                           exposure::Float64=1.0,
                           L_avg::Float64=1.0,
                           key::Float64=0.18)
    X, Y, Z = spectrum_to_XYZ(ν_grid, I_ν, cmfs)

    # Scale by exposure
    X *= exposure
    Y *= exposure
    Z *= exposure

    # Tone map the luminance (Y channel)
    Y_mapped = tone_map_reinhard(Y; key=key, L_avg=L_avg)
    scale = Y > 0 ? Y_mapped / Y : 0.0
    X_tm = X * scale
    Y_tm = Y * scale
    Z_tm = Z * scale

    # XYZ → linear sRGB
    R_lin, G_lin, B_lin = XYZ_to_linear_sRGB(X_tm, Y_tm, Z_tm)

    # Clamp negatives (out-of-gamut)
    R_lin = max(0.0, R_lin)
    G_lin = max(0.0, G_lin)
    B_lin = max(0.0, B_lin)

    # Gamma correction
    R_srgb = linear_sRGB_to_sRGB(R_lin)
    G_srgb = linear_sRGB_to_sRGB(G_lin)
    B_srgb = linear_sRGB_to_sRGB(B_lin)

    # To 8-bit
    r = clamp(round(Int, R_srgb * 255), 0, 255)
    g = clamp(round(Int, G_srgb * 255), 0, 255)
    b = clamp(round(Int, B_srgb * 255), 0, 255)

    return r, g, b
end

end # module
