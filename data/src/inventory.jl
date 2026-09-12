# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Record what each reanalysis netCDF file actually says, before anything trusts it: the
# variables, their unit strings verbatim, the time axis, the lattice, and the packing
# attributes. One row per data variable per file, written to `inventory.csv`.
#
#     julia --project data/src/inventory.jl [dir ...]
#
# With no arguments it inventories every `.nc` under `$AFRICA_DATADIR/era5`. The table is what a
# catalogue row is checked against, so it is the first thing to read when a file arrives.

module Inventory

using CSV
using DataFrames
using Dates
using NCDatasets
using Statistics

include("common.jl")

# Names a file may use for each axis; the first present wins.
const TIME_NAMES = ("valid_time", "time", "date")
const LAT_NAMES = ("latitude", "lat", "y")
const LON_NAMES = ("longitude", "lon", "x")

# The first of `names` that is a variable in `ds`, or `nothing`.
function _firstpresent(ds, names)
    for n in names
        haskey(ds, n) && return n
    end
    return nothing
end

# Spacing statistics of a coordinate vector: the step, and whether every step matches it to
# one part in a thousand. Gaussian latitudes fail this; regular lattices pass.
function _spacing(v::AbstractVector)
    length(v) < 2 && return (step = missing, regular = missing)
    d = diff(v)
    step = median(d)
    return (step = step, regular = all(abs.(d .- step) .<= abs(step) / 1000))
end

# The time axis as `DateTime`s, whatever calendar the file declares; `nothing` if there is none.
function _times(ds, name)
    isnothing(name) && return nothing
    raw = ds[name][:]
    return [DateTime(Dates.year(t), Dates.month(t), Dates.day(t)) for t in raw]
end

# A CF attribute as a string, or `missing`.
_attr(var, key) = haskey(var.attrib, key) ? string(var.attrib[key]) : missing

"""
    inventory_file(path)

Return a `DataFrame` with one row per data variable in the netCDF file at `path`: the file,
the variable, its dimensions, `units` verbatim, storage type, `scale_factor`, `add_offset`,
`_FillValue`, `grid_mapping`; the time axis calendar, first and last dates, number of slices,
spacing in days and whether it is even; and the latitude and longitude ranges, spacings and
whether each lattice is regular.
"""
function inventory_file(path::AbstractString)
    rows = DataFrame()
    NCDataset(path) do ds
        tname = _firstpresent(ds, TIME_NAMES)
        latname = _firstpresent(ds, LAT_NAMES)
        lonname = _firstpresent(ds, LON_NAMES)
        times = _times(ds, tname)
        tspace = isnothing(times) ? (step = missing, regular = missing) :
                 _spacing(Dates.value.(times) ./ 86_400_000)
        lat = isnothing(latname) ? nothing : Float64.(ds[latname][:])
        lon = isnothing(lonname) ? nothing : Float64.(ds[lonname][:])
        latspace = isnothing(lat) ? (step = missing, regular = missing) :
                   _spacing(lat)
        lonspace = isnothing(lon) ? (step = missing, regular = missing) :
                   _spacing(lon)
        axes = filter(!isnothing, (tname, latname, lonname))
        for name in keys(ds)
            name in axes && continue
            var = ds[name]
            # Coordinate-like helpers (bounds, numbers) carry no data of their own.
            ndims(var) < 2 && continue
            push!(rows,
                  (file = basename(path), variable = name,
                   dims = join(NCDatasets.dimnames(var), ","),
                   units = _attr(var, "units"),
                   storage = string(eltype(var.var)),
                   scale_factor = _attr(var, "scale_factor"),
                   add_offset = _attr(var, "add_offset"),
                   fillvalue = _attr(var, "_FillValue"),
                   grid_mapping = _attr(var, "grid_mapping"),
                   time_name = something(tname, missing),
                   calendar = isnothing(tname) ? missing :
                              _attr(ds[tname], "calendar"),
                   time_first = isnothing(times) ? missing : first(times),
                   time_last = isnothing(times) ? missing : last(times),
                   time_count = isnothing(times) ? missing : length(times),
                   time_step_days = tspace.step,
                   time_even = tspace.regular,
                   lat_min = isnothing(lat) ? missing : minimum(lat),
                   lat_max = isnothing(lat) ? missing : maximum(lat),
                   lat_step = latspace.step,
                   lat_regular = latspace.regular,
                   lon_min = isnothing(lon) ? missing : minimum(lon),
                   lon_max = isnothing(lon) ? missing : maximum(lon),
                   lon_step = lonspace.step,
                   lon_regular = lonspace.regular),
                  promote = true)
        end
    end
    return rows
end

"""
    inventory(dirs) -> DataFrame

Inventory every `.nc` file under each of `dirs`, recursively, one table for all of them.
"""
function inventory(dirs)
    tables = DataFrame[]
    for dir in dirs
        for (root, _, files) in walkdir(dir), f in files
            endswith(f, ".nc") || continue
            push!(tables, inventory_file(joinpath(root, f)))
        end
    end
    return isempty(tables) ? DataFrame() : reduce(vcat, tables, cols = :union)
end

# Write the inventory of the given directories (default `era5/` under the data directory) to
# `inventory.csv`, and say how many files it saw.
function main(argv)
    dirs = isempty(filter(a -> !startswith(a, "--"), argv)) ?
           [joinpath(datadir(), "era5")] :
           filter(a -> !startswith(a, "--"), argv)
    table = inventory(dirs)
    out = joinpath(datadir(), "inventory.csv")
    CSV.write(out, table)
    write_provenance(out,
                     source = "netCDF attributes of " *
                              join(repopath.(dirs), ", "),
                     files = isempty(table) ? String[] : unique(table.file))
    println("Inventoried ", isempty(table) ? 0 : length(unique(table.file)),
            " file(s), ", nrow(table), " variable(s) -> ", out)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Inventory.main(ARGS)
end
