# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Turn GBIF's simple CSV into the record table the species step and the burn-in use: one row per
# usable record, seven columns, streamed so the multi-gigabyte input never sits in memory.
#
#     julia --project data/src/gbif_clean.jl [era5|early|input.csv] [--force]
#
# The `era5` download (1940 onward) cleans to `records.csv`, which the milestone-1 pipeline
# reads; the `early` download (1806-1939) cleans to `records_early.csv`. Both when neither is
# named and both exist.
# Keeps a record when it has a species-level name and key, coordinates, a month and a year, an
# accepted basis of record, a coordinate uncertainty at or below `MAX_UNCERTAINTY_M` (or none
# stated), and does not share a 0.02 degree cell with a botanic garden or a country centroid,
# when `data/gardens.csv` and `data/centroids.csv` exist to say where those are. Writes
# `records.csv` and `records_by_decade.csv`, the latter for choosing a start year.

module GBIFClean

using CSV
using DataFrames

include("common.jl")

const MAX_UNCERTAINTY_M = 10_000.0
const CELL_DEG = 0.02
const KEEP_BASIS = Set(["PRESERVED_SPECIMEN", "HUMAN_OBSERVATION",
                           "OBSERVATION",
                           "MATERIAL_SAMPLE", "MATERIAL_CITATION", "OCCURRENCE",
                           "MACHINE_OBSERVATION"])
const COLUMNS = [:gbifID, :speciesKey, :species, :decimalLatitude,
    :decimalLongitude,
    :year, :month]

"""
    cellkey(lat, lon; step = CELL_DEG)

Return the index of the `step`-degree cell containing `(lat, lon)`, as a pair of integers from
the south-west corner of the globe. Two points share a cell when their keys are equal; this is
the anti-join key the garden and centroid filter uses.
"""
function cellkey(lat, lon; step = CELL_DEG)
    return (floor(Int, (lat + 90) / step),
            floor(Int, (lon + 180) / step))
end

# Latitude and longitude columns of a table of places, whatever they are called.
function _latlon(df::DataFrame)
    lower = Dict(lowercase(string(n)) => n for n in names(df))
    lat = get(lower, "latitude",
              get(lower, "decimallatitude",
                  get(lower, "lat", get(lower, "y", nothing))))
    lon = get(lower, "longitude",
              get(lower, "decimallongitude",
                  get(lower, "lon",
                      get(lower, "long", get(lower, "x", nothing)))))
    (isnothing(lat) || isnothing(lon)) &&
        error("cannot find latitude/longitude columns among " *
              join(names(df), ", "))
    return df[!, lat], df[!, lon]
end

"""
    excluded_cells(paths...) -> Set

Return the set of 0.02 degree cell keys occupied by the places listed in the CSV files at
`paths` (botanic gardens, country centroids); files that do not exist are skipped with a note.
"""
function excluded_cells(paths...)
    cells = Set{Tuple{Int, Int}}()
    for p in paths
        if !isfile(p)
            println("  no ", basename(p), ": that filter is skipped")
            continue
        end
        lat, lon = _latlon(CSV.read(p, DataFrame))
        for (a, b) in zip(lat, lon)
            (ismissing(a) || ismissing(b)) && continue
            push!(cells, cellkey(Float64(a), Float64(b)))
        end
    end
    return cells
end

# Whether one parsed GBIF row passes every filter. `row` is a `CSV.Row` or named tuple with the
# simple-CSV column names.
function keep(row, excluded; maxuncertainty = MAX_UNCERTAINTY_M)
    sp = row.species
    (ismissing(sp) || isempty(sp) || ismissing(row.speciesKey)) && return false
    (ismissing(row.decimalLatitude) || ismissing(row.decimalLongitude)) &&
        return false
    (ismissing(row.year) || ismissing(row.month)) && return false
    (1 <= row.month <= 12) || return false
    basis = row.basisOfRecord
    (ismissing(basis) || !(basis in KEEP_BASIS)) && return false
    u = row.coordinateUncertaintyInMeters
    (!ismissing(u) && u > maxuncertainty) && return false
    isempty(excluded) && return true
    return !(cellkey(Float64(row.decimalLatitude),
                     Float64(row.decimalLongitude)) in excluded)
end

"""
    clean_records(input, output; excluded = Set(), maxuncertainty = MAX_UNCERTAINTY_M)

Stream the GBIF simple CSV at `input` (tab-separated) and write the records that pass
[`keep`](@ref) to `output` with the columns `gbifID, speciesKey, species, decimalLatitude,
decimalLongitude, year, month`. Return `(kept, seen, bydecade)`, the last a `Dict` of decade
start year to records kept.
"""
function clean_records(input::AbstractString, output::AbstractString;
                       excluded = Set{Tuple{Int, Int}}(),
                       maxuncertainty = MAX_UNCERTAINTY_M)
    rows = CSV.Rows(input, delim = '\t', quoted = false,
                    select = [:gbifID, :speciesKey, :species, :decimalLatitude,
                        :decimalLongitude, :year, :month, :basisOfRecord,
                        :coordinateUncertaintyInMeters],
                    types = Dict(:gbifID => Int64, :speciesKey => Int64,
                                 :species => String,
                                 :decimalLatitude => Float64,
                                 :decimalLongitude => Float64, :year => Int,
                                 :month => Int, :basisOfRecord => String,
                                 :coordinateUncertaintyInMeters => Float64),
                    reusebuffer = true, silencewarnings = true)
    seen = 0
    kept = 0
    bydecade = Dict{Int, Int}()
    buffer = DataFrame([c => Vector{Any}() for c in COLUMNS])
    open(output, "w") do io
        CSV.write(io, buffer)   # the header
        for row in rows
            seen += 1
            keep(row, excluded, maxuncertainty = maxuncertainty) || continue
            kept += 1
            decade = 10 * fld(row.year, 10)
            bydecade[decade] = get(bydecade, decade, 0) + 1
            push!(buffer,
                  (row.gbifID, row.speciesKey, row.species, row.decimalLatitude,
                   row.decimalLongitude, row.year, row.month))
            if nrow(buffer) >= 100_000
                CSV.write(io, buffer, append = true)
                empty!(buffer)
            end
        end
        return CSV.write(io, buffer, append = true)
    end
    return (kept = kept, seen = seen, bydecade = bydecade)
end

# The downloads to clean as `(input, output)` pairs: a CSV named on the command line (to
# `records.csv`), or each of the named downloads present under `gbif/<name>/`, newest CSV first.
function _inputs(argv)
    named = filter(a -> !startswith(a, "--"), argv)
    out(name) = joinpath(datadir(),
                         name == "era5" ? "records.csv" : "records_$name.csv")
    pairs = Tuple{String, String}[]
    wanted = isempty(named) ? ["era5", "early"] : named
    for w in wanted
        if endswith(w, ".csv")
            push!(pairs, (w, out("era5")))
            continue
        end
        dir = joinpath(datadir(), "gbif", w)
        files = isdir(dir) ?
                filter(f -> endswith(f, ".csv"), readdir(dir, join = true)) :
                String[]
        if isempty(files)
            isempty(named) ||
                error("no GBIF CSV under $dir; run gbif_download.jl $w first")
            continue
        end
        push!(pairs, (last(sort(files, by = mtime)), out(w)))
    end
    isempty(pairs) &&
        error("no GBIF downloads under " * joinpath(datadir(), "gbif") *
              "; run gbif_download.jl first")
    return pairs
end

function main(argv)
    excluded = nothing
    for (input, output) in _inputs(argv)
        if isfile(output) && !forced(argv)
            println("Already have ", output,
                    "; nothing to do (use --force to rebuild).")
            continue
        end
        if isnothing(excluded)
            excluded = excluded_cells(joinpath(datadir(), "gardens.csv"),
                                      joinpath(datadir(), "centroids.csv"))
        end
        println("Cleaning ", input, " (", length(excluded),
                " excluded 0.02 degree cells) ...")
        result = clean_records(input, output, excluded = excluded)
        decades = sort(collect(result.bydecade))
        CSV.write(replace(output, ".csv" => "_by_decade.csv"),
                  DataFrame(decade = first.(decades), records = last.(decades)))
        write_provenance(output, source = "GBIF simple CSV " * input,
                         records_seen = result.seen, records_kept = result.kept,
                         excluded_cells = length(excluded),
                         max_uncertainty_m = MAX_UNCERTAINTY_M,
                         basis_of_record = sort(collect(KEEP_BASIS)))
        println("Kept ", result.kept, " of ", result.seen, " records -> ",
                output)
        for (d, n) in decades
            println("  ", d, "s: ", n)
        end
    end
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    GBIFClean.main(ARGS)
end
