# SPDX-License-Identifier: LGPL-3.0-or-later
#
# From cleaned GBIF records and the ERA5 monthly series, the per-species table the run is
# parameterised from, and the Africa subset of the records the burn-in seeds from.
#
#     julia --project data/src/species_traits.jl [--force]
#
# For every record: the ERA5 cell containing it and the twelve months ending in the month it was
# recorded; over that window the minimum monthly 2 m temperature, and the means of total
# precipitation, net solar radiation and soil water, as ClimatePref did for the published
# envelopes. Per species: the mean and standard deviation of those across its records within
# Africa. Records are inside Africa when they fall in an active cell of a geographic study grid
# cut to Natural Earth's continent outline, which needs no projection and follows the original's
# continent raster.
#
# Reads the netCDF files directly with NCDatasets rather than through EcoSISTEM's reader, because
# the join is at the data's own resolution and what the files call their axes is not yet known
# (see `inventory.jl`). Writes `species.csv`, `records_africa.csv` and `record_climate.csv`.
# Placing records on the run's grid and the effort raster are done at run time, since the grid
# is chosen then.

module SpeciesTraits

using Africa
using CSV
using DataFrames
using Dates
using EcoSISTEM
using EcoSISTEM.Units
using NCDatasets
using RasterDataSources
using Statistics
using Unitful
using Unitful.DefaultSymbols

include("common.jl")

const VARIABLES = ("t2m", "tp", "ssr", "swvl1")
const WINDOW = 12
const MIN_RECORDS = 2

# --- the climate stacks -------------------------------------------------------------------

"""
    ClimateStack

One variable's monthly series over the box: `lat` and `lon` coordinates (whatever order the
file stores them in), `months` as `(year, month)` pairs in time order, and `data` indexed
`[lon, lat, time]` as netCDF stores it, with missing values as `NaN`.
"""
struct ClimateStack
    variable::String
    lat::Vector{Float64}
    lon::Vector{Float64}
    months::Vector{Tuple{Int, Int}}
    data::Array{Float32, 3}
end

# Names a file may use for each axis; the first present wins (as `inventory.jl`).
const TIME_NAMES = ("valid_time", "time", "date")
const LAT_NAMES = ("latitude", "lat", "y")
const LON_NAMES = ("longitude", "lon", "x")

function _firstpresent(ds, names)
    for n in names
        haskey(ds, n) && return n
    end
    return error("none of " * join(names, ", ") * " is in the file")
end

# Read one decade file of `variable`: coordinates, months, and the data with the file's
# packing undone by NCDatasets and its fill turned into `NaN`.
function _readfile(path::AbstractString, variable::AbstractString)
    NCDataset(path) do ds
        lat = Float64.(ds[_firstpresent(ds, LAT_NAMES)][:])
        lon = Float64.(ds[_firstpresent(ds, LON_NAMES)][:])
        times = ds[_firstpresent(ds, TIME_NAMES)][:]
        months = [(Dates.year(t), Dates.month(t)) for t in times]
        raw = ds[variable][:, :, :]
        data = Array{Float32, 3}(undef, size(raw))
        @inbounds for i in eachindex(raw)
            v = raw[i]
            data[i] = ismissing(v) ? NaN32 : Float32(v)
        end
        return (lat = lat, lon = lon, months = months, data = data)
    end
end

"""
    load_stack(dir, variable) -> ClimateStack

Read every `era5_<variable>_*.nc` under `dir` and concatenate them in time order into one
[`ClimateStack`](@ref), checking that the files share one lattice and that the months run on
without a gap.
"""
function load_stack(dir::AbstractString, variable::AbstractString)
    files = sort(filter(f -> startswith(basename(f), "era5_$(variable)_") &&
                             endswith(f, ".nc"), readdir(dir, join = true)))
    isempty(files) && error("no era5_$(variable)_*.nc under $dir")
    parts = [_readfile(f, variable) for f in files]
    lat, lon = parts[1].lat, parts[1].lon
    for p in parts[2:end]
        (p.lat == lat && p.lon == lon) ||
            error("the $variable files are not on one lattice")
    end
    months = reduce(vcat, (p.months for p in parts))
    for i in 2:length(months)
        y, m = months[i - 1]
        months[i] == (m == 12 ? (y + 1, 1) : (y, m + 1)) ||
            error("the $variable series has a gap or repeat at $(months[i])")
    end
    data = cat((p.data for p in parts)..., dims = 3)
    return ClimateStack(variable, lat, lon, months, data)
end

# Index of the coordinate nearest `v` in `coords`, which may run either way.
function _nearest(coords::Vector{Float64}, v::Real)
    best = 1
    bestd = abs(coords[1] - v)
    @inbounds for i in 2:length(coords)
        d = abs(coords[i] - v)
        if d < bestd
            best = i
            bestd = d
        end
    end
    return best
end

# The time index of `(year, month)` in `stack`, or `nothing` if the series does not hold it.
function _timeindex(stack::ClimateStack, year::Integer, month::Integer)
    y0, m0 = stack.months[1]
    t = (year - y0) * 12 + (month - m0) + 1
    return 1 <= t <= length(stack.months) ? t : nothing
end

# --- per-record summaries ---------------------------------------------------------------

"""
    window_summary(stack, lat, lon, year, month; window = WINDOW, reducer = mean)

Return `reducer` over the `window` monthly values of `stack` ending in `(year, month)` at the
cell nearest `(lat, lon)`, ignoring `NaN`s; `NaN` when the window is not wholly inside the
series or has no data.
"""
function window_summary(stack::ClimateStack, lat::Real, lon::Real,
                        year::Integer,
                        month::Integer; window::Integer = WINDOW,
                        reducer = mean)
    t = _timeindex(stack, year, month)
    (isnothing(t) || t < window) && return NaN
    i = _nearest(stack.lon, lon)
    j = _nearest(stack.lat, lat)
    vals = @view stack.data[i, j, (t - window + 1):t]
    good = filter(!isnan, vals)
    return isempty(good) ? NaN : Float64(reducer(good))
end

"""
    record_climate(records, stacks) -> DataFrame

Return `records` with four columns added, the twelve-month summaries at each record's cell and
date: `tmin` (minimum monthly 2 m temperature, K), `tp` (mean total precipitation), `ssr` (mean
net solar radiation) and `swvl1` (mean volumetric soil water). `stacks` is a `Dict` from
variable name to [`ClimateStack`](@ref). Records whose window is not covered get `NaN`.
"""
function record_climate(records::DataFrame, stacks::Dict{String, ClimateStack})
    n = nrow(records)
    tmin = fill(NaN, n)
    tp = fill(NaN, n)
    ssr = fill(NaN, n)
    swvl1 = fill(NaN, n)
    lat, lon = records.decimalLatitude, records.decimalLongitude
    yr, mo = records.year, records.month
    Threads.@threads for r in 1:n
        tmin[r] = window_summary(stacks["t2m"], lat[r], lon[r], yr[r], mo[r],
                                 reducer = minimum)
        tp[r] = window_summary(stacks["tp"], lat[r], lon[r], yr[r], mo[r])
        ssr[r] = window_summary(stacks["ssr"], lat[r], lon[r], yr[r], mo[r])
        swvl1[r] = window_summary(stacks["swvl1"], lat[r], lon[r], yr[r], mo[r])
    end
    out = copy(records)
    out.tmin = tmin
    out.tp = tp
    out.ssr = ssr
    out.swvl1 = swvl1
    return out
end

# --- per-species table -----------------------------------------------------------------

# Standard deviation that is `NaN` rather than an error for fewer than two values.
_sd(v) = length(v) < 2 ? NaN : std(v)

"""
    species_table(climate; minrecords = MIN_RECORDS) -> DataFrame

Group the output of [`record_climate`](@ref) by species and return one row per species with
at least `minrecords` records carrying climate: `speciesKey, species, numrecords, tmin_mean,
tmin_sd, tp_mean, tp_sd, ssr_mean, swvl1_mean`. Means and standard deviations are over the
records with a value for that variable.
"""
function species_table(climate::DataFrame; minrecords::Integer = MIN_RECORDS)
    good = filter(r -> !isnan(r.tmin), climate)
    grouped = groupby(good, [:speciesKey, :species])
    table = combine(grouped,
                    nrow => :numrecords,
                    :tmin => mean => :tmin_mean,
                    :tmin => _sd => :tmin_sd,
                    :tp => (v -> mean(filter(!isnan, v))) => :tp_mean,
                    :tp => (v -> _sd(filter(!isnan, v))) => :tp_sd,
                    :ssr => (v -> mean(filter(!isnan, v))) => :ssr_mean,
                    :swvl1 => (v -> mean(filter(!isnan, v))) => :swvl1_mean)
    return sort!(filter(r -> r.numrecords >= minrecords, table), :species)
end

# --- the Africa subset -----------------------------------------------------------------

"""
    in_africa(records) -> BitVector

Return, for each record, whether it lies in an active cell of a geographic study grid cut to
Natural Earth's continent-level Africa: WorldClim's 10 arcminute annual mean temperature
positions the grid (cached after the first download) and the outline decides the cells. A
record over sea or outside the outline is `false`.
"""
function in_africa(records::DataFrame)
    temperature = SourceSpec(WorldClim{BioClim}, :bio1)
    # No `crs`: the grid stays WorldClim's own geographic lattice, which can be built and
    # indexed by latitude and longitude though not simulated. 10 arcminutes is WorldClim's
    # own resolution, so nothing is resampled. `Unitful.d` rather than `day`, which `Dates`
    # (loaded by `common.jl`) also exports.
    area = StudyArea(regime = temperature,
                     within = africa_within(outline = true),
                     cellsize = 10.0arcminute, verbosity = :silent)
    environment = GridHabitat(regime = temperature,
                              supply = UniformSpec(1.0kJ / (m^2 * Unitful.d),
                                                   axis = SolarRadiation),
                              area = area)
    active = parent(environment.active)
    ny, nx = size(active)
    inside = falses(nrow(records))
    for (r, (lat, lon)) in enumerate(zip(records.decimalLatitude,
                                         records.decimalLongitude))
        # `getcellat` is public, not exported. It throws for a place off the grid, which is
        # the one error expected here; anything else is a bug and must not be swallowed.
        cell = try
            EcoSISTEM.getcellat(environment, LatLong(lat * °, lon * °))
        catch e
            e isa UndefVarError && rethrow()
            nothing
        end
        isnothing(cell) && continue
        # The record's own cell, or any of its eight neighbours: a coastal record whose 10
        # arcminute cell is mostly sea sits in an inactive cell (Lagos does), and a one-cell
        # tolerance keeps it without reaching the offshore islands the outline excludes.
        j, i = Tuple(cell)
        inside[r] = any(active[jj, ii]
                        for jj in max(1, j - 1):min(ny, j + 1),
                            ii in max(1, i - 1):min(nx, i + 1))
    end
    return inside
end

# A `(year, month)` pair as `YYYY-MM`.
_monthlabel((y, m)) = string(y, "-", lpad(m, 2, "0"))

# --- the script ----------------------------------------------------------------------

function main(argv)
    out = joinpath(datadir(), "species.csv")
    if isfile(out) && !forced(argv)
        println("Already have ", out,
                "; nothing to do (use --force to rebuild).")
        return nothing
    end
    records = CSV.read(joinpath(datadir(), "records.csv"), DataFrame)
    println("Read ", nrow(records), " records")

    inside = in_africa(records)
    africa = records[inside, :]
    println(nrow(africa), " records inside Africa's outline")
    africapath = joinpath(datadir(), "records_africa.csv")
    CSV.write(africapath, africa)
    write_provenance(africapath,
                     source = "records.csv within Natural Earth Africa",
                     records = nrow(africa))

    era5 = joinpath(datadir(), "era5")
    stacks = Dict(v => load_stack(era5, v) for v in VARIABLES)
    first_month, last_month = stacks["t2m"].months[1], stacks["t2m"].months[end]
    println("ERA5 months ", first_month, " to ", last_month)

    climate = record_climate(africa, stacks)
    climatepath = joinpath(datadir(), "record_climate.csv")
    CSV.write(climatepath, climate)
    write_provenance(climatepath, source = "records_africa.csv joined to era5/",
                     window_months = WINDOW,
                     records_with_climate = count(!isnan, climate.tmin))

    table = species_table(climate)
    CSV.write(out, table)
    write_provenance(out, source = "record_climate.csv grouped by species",
                     species = nrow(table), min_records = MIN_RECORDS,
                     era5_first_month = _monthlabel(first_month),
                     era5_last_month = _monthlabel(last_month))
    println(nrow(table), " species with at least ", MIN_RECORDS,
            " records carrying climate -> ", out)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    SpeciesTraits.main(ARGS)
end
