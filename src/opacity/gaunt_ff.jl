#=
Free-free Gaunt factor interpolation from Sutherland (1998) table.

Source: Sutherland, MNRAS 300, 321 (1998).
Table: gffgu.dat distributed with McPHAC (refs/code/McPHAC/gffgu.dat).

File format: 3 columns per row — γ², u, g̃_ff
where u = hν/(k_BT), γ² = Z²Ry/(k_BT), Ry = 13.6057 eV.
Note: McPHAC GauntFactor.c reads col1 as g2str and col2 as ustr,
then stores gauntff[i][0]=u, gauntff[i][1]=g2 (lines 37-38).
=#

module GauntFactor

using ..PhysicalConstants: h, k_B, eV

export load_gaunt_table, gaunt_ff, GauntTable

"Rydberg energy [erg]. CODATA 2018."
const Ry_erg = 13.605693122994 * eV

"""
Structured Gaunt factor table for bilinear interpolation.
"""
struct GauntTable
    u_grid::Vector{Float64}    # sorted unique u values
    g2_grid::Vector{Float64}   # sorted unique γ² values
    gff::Matrix{Float64}       # gff[iu, ig2]
end

"""
    load_gaunt_table(path) → GauntTable

Load the Sutherland (1998) gffgu.dat table.
"""
function load_gaunt_table(path::String)::GauntTable
    g2_raw = Float64[]
    u_raw = Float64[]
    gff_raw = Float64[]

    for line in eachline(path)
        parts = split(strip(line))
        length(parts) >= 3 || continue
        push!(g2_raw, parse(Float64, parts[1]))   # col 1 = γ²
        push!(u_raw, parse(Float64, parts[2]))     # col 2 = u
        push!(gff_raw, parse(Float64, parts[3]))   # col 3 = g̃_ff
    end

    @assert length(g2_raw) > 100 "Too few Gaunt data: $(length(g2_raw))"

    # Extract unique sorted grids
    u_grid = sort(unique(u_raw))
    g2_grid = sort(unique(g2_raw))
    Nu = length(u_grid)
    Ng = length(g2_grid)

    # Build 2D array
    gff = zeros(Nu, Ng)
    for i in eachindex(g2_raw)
        iu = searchsortedfirst(u_grid, u_raw[i])
        ig = searchsortedfirst(g2_grid, g2_raw[i])
        if iu <= Nu && ig <= Ng
            gff[iu, ig] = gff_raw[i]
        end
    end

    @assert all(gff .> 0) "Gaunt table has zero/negative entries"
    return GauntTable(u_grid, g2_grid, gff)
end

"""
    gaunt_ff(ν, T, table) → g̃_ff

Thermally-averaged free-free Gaunt factor via bilinear interpolation.

Arguments:
- `ν`: photon frequency [Hz]
- `T`: temperature [K]
- `table`: loaded GauntTable
"""
function gaunt_ff(ν::Float64, T::Float64, table::GauntTable)::Float64
    @assert ν > 0 && T > 0

    u = h * ν / (k_B * T)
    g2 = Ry_erg / (k_B * T)   # Z=1 for hydrogen

    return _interp_bilinear(u, g2, table)
end

"""
Bilinear interpolation with boundary clamping (matches McPHAC behaviour).
"""
function _interp_bilinear(u::Float64, g2::Float64, t::GauntTable)::Float64
    ug = t.u_grid
    gg = t.g2_grid
    Nu = length(ug)
    Ng = length(gg)

    # Clamp to table range (McPHAC boundary behaviour)
    u_c = clamp(u, ug[1], ug[end])
    g2_c = clamp(g2, gg[1], gg[end])

    # Find bracketing indices
    iu = _bracket(ug, u_c)
    ig = _bracket(gg, g2_c)

    # Interpolation fractions
    s = iu < Nu ? (u_c - ug[iu]) / (ug[iu+1] - ug[iu]) : 0.0
    r = ig < Ng ? (g2_c - gg[ig]) / (gg[ig+1] - gg[ig]) : 0.0

    iu2 = min(iu + 1, Nu)
    ig2 = min(ig + 1, Ng)

    # Bilinear
    gff = (1-s)*(1-r)*t.gff[iu, ig] +
          s*(1-r)*t.gff[iu2, ig] +
          s*r*t.gff[iu2, ig2] +
          (1-s)*r*t.gff[iu, ig2]

    @assert gff > 0 "Gaunt factor non-positive: $gff at u=$u, g2=$g2"
    return gff
end

"Find index i such that grid[i] ≤ x < grid[i+1], clamped to [1, N-1]."
function _bracket(grid::Vector{Float64}, x::Float64)::Int
    i = searchsortedlast(grid, x)
    return clamp(i, 1, length(grid) - 1)
end

end # module
