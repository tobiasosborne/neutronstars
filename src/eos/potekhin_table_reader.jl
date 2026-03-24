#=
Parser for Potekhin magnetised hydrogen EOS + opacity tables.

Reads the .dat files from refs/potekhin_tables/ which contain 14 columns:
  EOS (cols 2-8), ionisation fractions (cols 9-12), Rosseland opacities (cols 13-14).

Table format documented in refs/potekhin_tables/hmagtab.txt.
Three file prefixes with different headers:
  hmm* (2-line header, lgB ≤ 12.0)
  hmag* (5-line header, 11.9 ≤ lgB ≤ 13.5)
  hmn* (5-line header, 13.5 ≤ lgB ≤ 15.0)
=#

module PotekhinTableReader

export PotekhinTableEntry, PotekhinIsotherm, PotekhinTable
export load_table, load_all_tables, merge_tables, lgB_from_filename

"""
Single row of a Potekhin table at fixed (lgT, lgB, lgR).
"""
struct PotekhinTableEntry
    lgR::Float64        # col 1:  lg(R), R = ρ/T₆³
    lgP::Float64        # col 2:  lg(P/bar)
    PV_NkT::Float64     # col 3:  dimensionless pressure
    U_NkT::Float64      # col 4:  dimensionless internal energy
    S_Nk::Float64       # col 5:  dimensionless entropy
    Cv_Nk::Float64      # col 6:  reduced heat capacity
    chit::Float64       # col 7:  (d ln P / d ln T)_V
    chir::Float64       # col 8:  (d ln P / d ln ρ)_T
    x_H::Float64        # col 9:  total atomic fraction
    x_H0::Float64       # col 10: ground-state atomic fraction
    x_H2::Float64       # col 11: molecular fraction
    x_pert::Float64     # col 12: perturbed fraction
    lgK0::Float64       # col 13: lg(κ_∥) [cm²/g]
    lgK1::Float64       # col 14: lg(κ_⊥) [cm²/g]
end

"""
One isotherm block: all density points at fixed (lgT, lgB).
"""
struct PotekhinIsotherm
    lgT::Float64
    lgB::Float64
    entries::Vector{PotekhinTableEntry}
end

"""
Full table from one .dat file (one B value, multiple T values).
"""
struct PotekhinTable
    lgB::Float64
    isotherms::Vector{PotekhinIsotherm}
end

"""
    lgB_from_filename(filename) → Float64

Extract lg(B) from a table filename.
  hmag12_5.dat → 12.5
  hmm10_5.dat  → 10.5
  hmm10a8.dat  → 10.845
  hmm11a8.dat  → 11.845
  hmn14_0.dat  → 14.0
"""
function lgB_from_filename(filename::AbstractString)::Float64
    base = basename(filename)
    base = replace(base, ".dat" => "")

    # Special cases: hmm10a8 → 10.845, hmm11a8 → 11.845
    if occursin("10a8", base)
        return 10.845
    elseif occursin("11a8", base)
        return 11.845
    end

    # Strip prefix (hmm, hmag, hmn)
    for prefix in ("hmag", "hmm", "hmn")
        if startswith(base, prefix)
            numstr = base[length(prefix)+1:end]
            # "12_5" → 12.5, "15_0" → 15.0
            numstr = replace(numstr, "_" => ".")
            return parse(Float64, numstr)
        end
    end

    error("Cannot parse lgB from filename: $filename")
end

"""
    load_table(filepath) → PotekhinTable

Parse a single .dat file into structured data.
"""
function load_table(filepath::AbstractString)::PotekhinTable
    lines = readlines(filepath)
    expected_lgB = lgB_from_filename(filepath)

    isotherms = PotekhinIsotherm[]
    current_lgT = NaN
    current_lgB = NaN
    current_entries = PotekhinTableEntry[]

    for line in lines
        stripped = strip(line)
        isempty(stripped) && continue

        tokens = split(stripped)
        ntok = length(tokens)

        # Try to parse as data row (14 numeric fields, first is lgR ∈ [-8, 4])
        if ntok >= 14
            vals = tryparse_floats(tokens[1:14])
            if vals !== nothing && -8.5 < vals[1] < 4.5
                entry = PotekhinTableEntry(
                    vals[1],  vals[2],  vals[3],  vals[4],
                    vals[5],  vals[6],  vals[7],  vals[8],
                    vals[9],  vals[10], vals[11], vals[12],
                    vals[13], vals[14]
                )
                push!(current_entries, entry)
                continue
            end
        end

        # Try to parse as lgT/lgB header (exactly 2 numeric fields)
        if ntok == 2
            vals = tryparse_floats(tokens)
            if vals !== nothing && 4.0 < vals[1] < 8.5 && 9.5 < vals[2] < 16.0
                # Save previous isotherm if any
                if !isempty(current_entries)
                    push!(isotherms, PotekhinIsotherm(current_lgT, current_lgB, current_entries))
                    current_entries = PotekhinTableEntry[]
                end
                current_lgT = vals[1]
                current_lgB = vals[2]
                continue
            end
        end

        # Otherwise: header/comment line, skip
    end

    # Save last isotherm
    if !isempty(current_entries)
        push!(isotherms, PotekhinIsotherm(current_lgT, current_lgB, current_entries))
    end

    @assert !isempty(isotherms) "No isotherms found in $filepath"

    # Verify lgB from header matches filename
    file_lgB = isotherms[1].lgB
    @assert abs(file_lgB - expected_lgB) < 0.01 "lgB mismatch in $filepath: filename=$expected_lgB, header=$file_lgB"

    # Verify all isotherms have the same lgB
    for iso in isotherms
        @assert abs(iso.lgB - file_lgB) < 0.01 "Inconsistent lgB in $filepath: $(iso.lgB) vs $file_lgB"
    end

    return PotekhinTable(file_lgB, isotherms)
end

"""
    load_all_tables(dir) → Vector{PotekhinTable}

Load all .dat files from a directory. Returns a vector (not a dict) because
some lgB values have two files with different T ranges:
  hmm12_0.dat (lgT 4.9-7.0, step 0.1) and hmag12_0.dat (lgT 5.3-7.0, step 0.05)
  hmag13_5.dat (lgT 5.3-7.0) and hmn13_5.dat (lgT 5.7-7.6)
"""
function load_all_tables(dir::AbstractString)::Vector{PotekhinTable}
    files = sort(filter(f -> endswith(f, ".dat"), readdir(dir)))
    @assert !isempty(files) "No .dat files found in $dir"

    tables = PotekhinTable[]
    for f in files
        path = joinpath(dir, f)
        table = load_table(path)
        push!(tables, table)
    end

    return tables
end

"""
    merge_tables(tables) → Dict{Float64, PotekhinTable}

Merge tables with the same lgB by combining their isotherms.
Where both tables have the same (lgT, lgB), keeps the one from the
file with more isotherms (finer T grid).
"""
function merge_tables(tables::Vector{PotekhinTable})::Dict{Float64, PotekhinTable}
    grouped = Dict{Float64, Vector{PotekhinTable}}()
    for t in tables
        key = round(t.lgB, digits=3)
        push!(get!(grouped, key, PotekhinTable[]), t)
    end

    merged = Dict{Float64, PotekhinTable}()
    for (lgB, group) in grouped
        if length(group) == 1
            merged[lgB] = group[1]
        else
            # Merge isotherms: collect all unique lgT values
            all_iso = Dict{Float64, PotekhinIsotherm}()
            # Process tables with fewer isotherms first so finer-grid ones overwrite
            sort!(group, by=t -> length(t.isotherms))
            for t in group
                for iso in t.isotherms
                    all_iso[round(iso.lgT, digits=3)] = iso
                end
            end
            sorted_iso = sort(collect(values(all_iso)), by=iso -> iso.lgT)
            merged[lgB] = PotekhinTable(lgB, sorted_iso)
        end
    end

    return merged
end

# --- Internal helpers ---

"""Try to parse a vector of string tokens as Float64s. Returns nothing on failure."""
function tryparse_floats(tokens)::Union{Nothing, Vector{Float64}}
    vals = Float64[]
    for t in tokens
        v = tryparse(Float64, t)
        v === nothing && return nothing
        push!(vals, v)
    end
    return vals
end

end # module
