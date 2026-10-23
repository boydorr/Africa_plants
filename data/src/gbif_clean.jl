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
#
# Keeps a record when it has a species-level name and key, coordinates, a month and a year, an
# accepted basis of record and a coordinate uncertainty at or below `MAX_UNCERTAINTY_M` (or none
# stated); when it is not within `PLACE_BUFFER_KM` of a country centroid listed in
# `data/centroids.csv` (which `centroids.jl` writes); and when its locality does not say it grew
# under glass. Three further properties are written as flag columns rather than applied as
# filters, since whether they disqualify a record depends on what it is used for: `on_lattice`
# (both coordinates on a quarter-degree lattice), `near_garden` (within `PLACE_BUFFER_KM` of a
# botanic garden listed in `data/gardens.csv`, which `gardens.jl` writes) and `garden_named` (the
# locality names a botanic garden, a nursery or cultivation). Writes the table and
# `<table>_by_decade.csv`, the latter for choosing a start year.

module GBIFClean

using CSV
using DataFrames

include("common.jl")

const MAX_UNCERTAINTY_M = 10_000.0
const PLACE_BUFFER_KM = 1.0
const LATTICE_DEG = 0.25
const KEEP_BASIS = Set(["PRESERVED_SPECIMEN", "HUMAN_OBSERVATION",
                           "OBSERVATION",
                           "MATERIAL_SAMPLE", "MATERIAL_CITATION", "OCCURRENCE",
                           "MACHINE_OBSERVATION"])
const COLUMNS = [:gbifID, :speciesKey, :species, :decimalLatitude,
    :decimalLongitude, :year, :month, :on_lattice, :near_garden, :garden_named]

# A locality saying the plant grew under glass, whose climate is not the site's.
const GREENHOUSE = r"greenhouse|glasshouse|hothouse|shade ?house"i

# A locality naming a botanic garden, a nursery or cultivation. "Arboretum" is left out, since in
# Africa it mostly names a forest reserve, and so is a bare "cultivated", mostly habitat text.
const GARDEN_NAMED = r"\bhort\b\.|\bhortus\b|botanic(al)?\s+gardens?|botan\.\s+garden|botaniese\s+tuin|botanischer\s+garten|jard[ií]n\s+bot|jardim\s+bot|\bnursery\b|\bkwekery\b|\bin cultivation\b"i

"""
    PlaceIndex(lat, lon; buffer_km = PLACE_BUFFER_KM)

A set of places, botanic gardens or country centroids, each excluding the records within
`buffer_km` of it. `lat` and `lon` are the places' coordinates in degrees. The places are held
in buckets of `step` degrees, so that a query compares a record only with the places near it.
"""
struct PlaceIndex
    lat::Vector{Float64}
    lon::Vector{Float64}
    buffer_km::Float64
    step::Float64
    buckets::Dict{Tuple{Int, Int}, Vector{Int}}
end

function PlaceIndex(lat::AbstractVector{<:Real}, lon::AbstractVector{<:Real};
                    buffer_km::Real = PLACE_BUFFER_KM)
    length(lat) == length(lon) ||
        error("latitude and longitude lists differ in length")
    step = 0.05
    buckets = Dict{Tuple{Int, Int}, Vector{Int}}()
    for (i, (a, b)) in enumerate(zip(lat, lon))
        push!(get!(buckets, (floor(Int, a / step), floor(Int, b / step)),
                   Int[]),
              i)
    end
    return PlaceIndex(Float64.(lat), Float64.(lon), Float64(buffer_km), step,
                      buckets)
end

Base.length(places::PlaceIndex) = length(places.lat)
Base.isempty(places::PlaceIndex) = isempty(places.lat)

function Base.show(io::IO, places::PlaceIndex)
    return print(io, "PlaceIndex(", length(places), " places, ",
                 places.buffer_km, " km)")
end

function Base.show(io::IO, ::MIME"text/plain", places::PlaceIndex)
    return print(io, "PlaceIndex of ", length(places),
                 " places excluding records within ", places.buffer_km, " km")
end

# == Functions ==

# Great-circle distance in km between two points given in degrees.
function _distance_km(lat1::Real, lon1::Real, lat2::Real, lon2::Real)
    h = sind((lat2 - lat1) / 2)^2 +
        cosd(lat1) * cosd(lat2) * sind((lon2 - lon1) / 2)^2
    return 2 * 6371.0 * asin(sqrt(min(h, 1.0)))
end

"""
    near(places, lat, lon)

Return whether `(lat, lon)` lies within the buffer of any of `places`. The buckets searched
cover the buffer in longitude at the point's own latitude, so a place just across a bucket edge
is still found.
"""
function near(places::PlaceIndex, lat::Real, lon::Real)
    isempty(places) && return false
    dlat = places.buffer_km / 110.57
    dlon = dlat / max(cosd(lat), 0.01)
    s = places.step
    for i in floor(Int, (lat - dlat) / s):floor(Int, (lat + dlat) / s),
        j in floor(Int, (lon - dlon) / s):floor(Int, (lon + dlon) / s)
        for k in get(places.buckets, (i, j), Int[])
            _distance_km(lat, lon, places.lat[k], places.lon[k]) <=
            places.buffer_km && return true
        end
    end
    return false
end

"""
    on_lattice(lat, lon; step = LATTICE_DEG)

Return whether both coordinates are whole multiples of `step` degrees. With the default quarter
degree these are whole and half degrees and the centres of one- and half-degree squares: where
a record lands when its position was set to a coarse grid cell or rounded rather than measured.
Centres of quarter-degree squares, such as `(-33.875, 18.875)`, are not multiples and pass.
"""
function on_lattice(lat::Real, lon::Real; step::Real = LATTICE_DEG)
    return _multiple(lat, step) && _multiple(lon, step)
end

# Whether `x` is a whole multiple of `step`, to within rounding in the decimal text GBIF writes.
_multiple(x::Real, step::Real) = abs(x / step - round(x / step)) < 1e-6

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
    load_places(paths...; buffer_km = PLACE_BUFFER_KM)

Return a [`PlaceIndex`](@ref) of the places listed in the CSV files at `paths` (country
centroids, botanic gardens), each marking records within `buffer_km`; files that do not exist are
skipped with a note.
"""
function load_places(paths::AbstractString...;
                     buffer_km::Real = PLACE_BUFFER_KM)
    lat = Float64[]
    lon = Float64[]
    for p in paths
        if !isfile(p)
            println("  no ", basename(p), ": that filter is skipped")
            continue
        end
        a, b = _latlon(CSV.read(p, DataFrame))
        for (x, y) in zip(a, b)
            (ismissing(x) || ismissing(y)) && continue
            push!(lat, x)
            push!(lon, y)
        end
    end
    return PlaceIndex(lat, lon, buffer_km = buffer_km)
end

# Whether a locality text, which may be missing, matches `pattern`.
function _mentions(locality, pattern::Regex)
    return !ismissing(locality) && occursin(pattern, locality)
end

"""
    greenhouse(locality)

Return whether a record's `locality` text says the plant grew under glass (`GREENHOUSE`); a
missing locality says nothing.
"""
greenhouse(locality) = _mentions(locality, GREENHOUSE)

"""
    garden_named(locality)

Return whether a record's `locality` text names a botanic garden, a nursery or cultivation
(`GARDEN_NAMED`); a missing locality names none.
"""
garden_named(locality) = _mentions(locality, GARDEN_NAMED)

# Whether one parsed GBIF row is usable before its position is judged: a species, coordinates, a
# month and year, an accepted basis of record and a small enough stated uncertainty. `row` is a
# `CSV.Row` or named tuple with the simple-CSV column names.
function _usable(row; maxuncertainty = MAX_UNCERTAINTY_M)
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
    return true
end

"""
    keep(row, centroids; maxuncertainty = MAX_UNCERTAINTY_M)

Return whether one parsed GBIF row is kept: usable, not near any of `centroids`, and not grown
under glass. `row` is a `CSV.Row` or named tuple with the simple-CSV column names, `locality`
among them.
"""
function keep(row, centroids::PlaceIndex; maxuncertainty = MAX_UNCERTAINTY_M)
    _usable(row, maxuncertainty = maxuncertainty) || return false
    near(centroids, row.decimalLatitude, row.decimalLongitude) && return false
    return !greenhouse(row.locality)
end

"""
    clean_records(input, output; centroids = PlaceIndex(Float64[], Float64[]),
                  gardens = PlaceIndex(Float64[], Float64[]),
                  maxuncertainty = MAX_UNCERTAINTY_M)

Stream the GBIF simple CSV at `input` (tab-separated) and write the records that pass
[`keep`](@ref) against `centroids` to `output` with the columns `COLUMNS`: the record, and the
`on_lattice`, `near_garden` (near any of `gardens`) and `garden_named` flags. Return `(kept, seen,
centroid, greenhouse, lattice, garden, named, bydecade)`: records kept and read, usable records
dropped near a centroid and under glass, kept records flagged on the lattice, near a garden and
naming one, and a `Dict` of decade start year to records kept.
"""
function clean_records(input::AbstractString, output::AbstractString;
                       centroids::PlaceIndex = PlaceIndex(Float64[], Float64[]),
                       gardens::PlaceIndex = PlaceIndex(Float64[], Float64[]),
                       maxuncertainty = MAX_UNCERTAINTY_M)
    rows = CSV.Rows(input, delim = '\t', quoted = false,
                    select = [:gbifID, :speciesKey, :species, :decimalLatitude,
                        :decimalLongitude, :year, :month, :basisOfRecord,
                        :coordinateUncertaintyInMeters, :locality],
                    types = Dict(:gbifID => Int64, :speciesKey => Int64,
                                 :species => String,
                                 :decimalLatitude => Float64,
                                 :decimalLongitude => Float64, :year => Int,
                                 :month => Int, :basisOfRecord => String,
                                 :coordinateUncertaintyInMeters => Float64,
                                 :locality => String),
                    reusebuffer = true, on_error = :collect)
    seen = 0
    kept = 0
    ncentroid = 0
    nglass = 0
    nlattice = 0
    ngarden = 0
    nnamed = 0
    bydecade = Dict{Int, Int}()
    buffer = DataFrame([c => Vector{Any}() for c in COLUMNS])
    open(output, "w") do io
        CSV.write(io, buffer)   # the header
        for row in rows
            seen += 1
            _usable(row, maxuncertainty = maxuncertainty) || continue
            lat, lon = row.decimalLatitude, row.decimalLongitude
            if near(centroids, lat, lon)
                ncentroid += 1
                continue
            end
            if greenhouse(row.locality)
                nglass += 1
                continue
            end
            lattice = on_lattice(lat, lon)
            garden = near(gardens, lat, lon)
            named = garden_named(row.locality)
            nlattice += lattice
            ngarden += garden
            nnamed += named
            kept += 1
            decade = 10 * fld(row.year, 10)
            bydecade[decade] = get(bydecade, decade, 0) + 1
            push!(buffer,
                  (row.gbifID, row.speciesKey, row.species, lat, lon, row.year,
                   row.month, lattice, garden, named))
            if nrow(buffer) >= 100_000
                CSV.write(io, buffer, append = true)
                empty!(buffer)
            end
        end
        return CSV.write(io, buffer, append = true)
    end
    return (kept = kept, seen = seen, centroid = ncentroid, greenhouse = nglass,
            lattice = nlattice, garden = ngarden, named = nnamed,
            bydecade = bydecade)
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
    centroidfile = joinpath(datadir(), "centroids.csv")
    gardenfile = joinpath(datadir(), "gardens.csv")
    centroids = nothing
    gardens = nothing
    for (input, output) in _inputs(argv)
        if isfile(output) && !forced(argv)
            println("Already have ", output,
                    "; nothing to do (use --force to rebuild).")
            continue
        end
        if isnothing(centroids)
            centroids = load_places(centroidfile)
            gardens = load_places(gardenfile)
        end
        println("Cleaning ", input, " (", length(centroids),
                " country centroids, ", length(gardens), " gardens) ...")
        result = clean_records(input, output, centroids = centroids,
                               gardens = gardens)
        decades = sort(collect(result.bydecade))
        CSV.write(replace(output, ".csv" => "_by_decade.csv"),
                  DataFrame(decade = first.(decades), records = last.(decades)))
        write_sidecar(output,
                      source = "GBIF simple CSV " * repopath(input),
                      records_seen = result.seen, records_kept = result.kept,
                      dropped_near_country_centroid = result.centroid,
                      dropped_greenhouse = result.greenhouse,
                      flagged_on_lattice = result.lattice,
                      flagged_near_garden = result.garden,
                      flagged_garden_named = result.named,
                      lattice_deg = LATTICE_DEG,
                      place_buffer_km = PLACE_BUFFER_KM,
                      place_files = [repopath(f)
                                     for f in (centroidfile, gardenfile)
                                     if isfile(f)],
                      max_uncertainty_m = MAX_UNCERTAINTY_M,
                      basis_of_record = sort(collect(KEEP_BASIS)))
        println("Kept ", result.kept, " of ", result.seen, " records -> ",
                output, " (dropped ", result.centroid,
                " near a country centroid, ",
                result.greenhouse, " under glass; flagged ", result.lattice,
                " on the lattice, ", result.garden, " near a garden, ",
                result.named, " naming one)")
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
