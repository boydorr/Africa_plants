# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Average the monthly MOD13C2 files `evi_download.jl` fetched into one EVI layer: at each 0.05
# degree cell, the mean of the monthly EVI values the product rates good or marginal, over every
# month from `EVIDownload.FIRST_MONTH` to `EVIDownload.LAST_MONTH`. The mean stands in for how many
# plants a place holds: collection effort in a cell is its record count over its EVI, so a record in
# a greener or less sampled cell seeds more individuals. Months the product rates snow or ice,
# cloudy or estimated are left out, and so are fill values.
#
#     julia --project data/src/evi_mean.jl [--force]
#
# Writes `$AFRICA_DATADIR/evi/evi_mean_2001_2025.nc`: global, `evi` (the mean, `NaN` where no month
# qualifies) and `months` (how many months went into it), on longitude and latitude cell centres.

module EVIMean

using ArchGDAL
using Dates
using NCDatasets

include("common.jl")
include("evi_download.jl")

const GRID = "MOD_Grid_monthly_CMG_VI"
const EVI_LAYER = "CMG 0.05 Deg Monthly EVI"
const RELIABILITY_LAYER = "CMG 0.05 Deg Monthly pixel reliability"
const EVI_FILL = Int16(-3000)
const EVI_SCALE = 10_000
const KEEP_RELIABILITY = (Int8(0), Int8(1))
const NX = 7200
const NY = 3600
const CELL = 0.05
const OUT = "evi_mean_2001_2025.nc"

"""
    accumulate!(total, count, evi, reliability)

Add one month to a running mean: at every cell whose `reliability` is good (0) or marginal (1) and
whose `evi` is not the fill value, add the EVI, divided by its storage scale, to `total` and one to
`count`. All four arrays have the same size.
"""
function accumulate!(total::AbstractArray{Float64}, count::AbstractArray{Int16},
                     evi::AbstractArray{<:Integer},
                     reliability::AbstractArray{<:Integer})
    size(total) == size(count) == size(evi) == size(reliability) ||
        error("the EVI, reliability and running-mean arrays differ in size")
    @inbounds for i in eachindex(evi)
        (reliability[i] in KEEP_RELIABILITY && evi[i] != EVI_FILL) || continue
        total[i] += evi[i] / EVI_SCALE
        count[i] += 1
    end
    return total, count
end

"""
    finish(total, count)

Return the mean EVI from a running `total` and `count` as `Float32`, `NaN` where no month counted.
"""
function finish(total::AbstractArray{Float64}, count::AbstractArray{Int16})
    return Float32[c == 0 ? NaN32 : Float32(t / c)
                   for (t, c) in zip(total, count)]
end

# The GDAL name of one layer of a MOD13C2 file.
_subdataset(path, layer) = "HDF4_EOS:EOS_GRID:\"$path\":$GRID:$layer"

# One band of a MOD13C2 layer, refusing a file not on the 0.05 degree global grid the output is
# written on. GDAL returns it `width x height`, the first row at the north.
function _readlayer(path::AbstractString, layer::AbstractString)
    return ArchGDAL.read(_subdataset(path, layer)) do ds
        gt = ArchGDAL.getgeotransform(ds)
        (ArchGDAL.width(ds), ArchGDAL.height(ds)) == (NX, NY) &&
        isapprox(gt, [-180.0, CELL, 0.0, 90.0, 0.0, -CELL]) ||
            error("$(basename(path)) is not on the 0.05 degree global grid: " *
                  "$(ArchGDAL.width(ds)) x $(ArchGDAL.height(ds)), geotransform $gt")
        return ArchGDAL.read(ArchGDAL.getband(ds, 1))
    end
end

"""
    mean_evi(paths)

Return `(evi, months)` over the MOD13C2 files at `paths`: the mean of each cell's good and
marginal monthly EVI as `Float32` (`NaN` where none), and how many months went into it, both
`NX x NY` with the first row at the north.
"""
function mean_evi(paths::AbstractVector{<:AbstractString})
    total = zeros(Float64, NX, NY)
    count = zeros(Int16, NX, NY)
    for (k, path) in enumerate(paths)
        accumulate!(total, count, _readlayer(path, EVI_LAYER),
                    _readlayer(path, RELIABILITY_LAYER))
        (k % 12 == 0 || k == length(paths)) &&
            (println("  ", k, " of ", length(paths), " months"); flush(stdout))
    end
    return (evi = finish(total, count), months = count)
end

"""
    write_mean(path, evi, months)

Write the mean EVI layer and its month counts to the netCDF file `path`, on longitude and latitude
cell centres of the 0.05 degree global grid.
"""
function write_mean(path::AbstractString, evi::AbstractMatrix{Float32},
                    months::AbstractMatrix{Int16})
    NCDataset(path, "c") do ds
        defDim(ds, "lon", NX)
        defDim(ds, "lat", NY)
        lon = defVar(ds, "lon", Float64, ("lon",),
                     attrib = ["units" => "degrees_east",
                         "standard_name" => "longitude"])
        lat = defVar(ds, "lat", Float64, ("lat",),
                     attrib = ["units" => "degrees_north",
                         "standard_name" => "latitude"])
        lon[:] = -180 .+ CELL .* ((1:NX) .- 0.5)
        lat[:] = 90 .- CELL .* ((1:NY) .- 0.5)
        v = defVar(ds, "evi", Float32, ("lon", "lat"), fillvalue = NaN32,
                   deflatelevel = 4,
                   attrib = ["units" => "1",
                       "long_name" => "mean enhanced vegetation index, good and marginal months"
                   ])
        v[:, :] = evi
        m = defVar(ds, "months", Int16, ("lon", "lat"), deflatelevel = 4,
                   attrib = ["units" => "1",
                       "long_name" => "number of months averaged"])
        m[:, :] = months
        return nothing
    end
    return path
end

function main(argv)
    dir = joinpath(datadir(), "evi")
    out = joinpath(dir, OUT)
    if isfile(out) && !forced(argv)
        println("Already have ", out,
                "; nothing to do (use --force to rebuild).")
        return nothing
    end
    names = filter(f -> startswith(f, "MOD13C2.A") && endswith(f, ".hdf"),
                   readdir(dir))
    first, last = EVIDownload.FIRST_MONTH, EVIDownload.LAST_MONTH
    names = filter(n -> first <= EVIDownload.granule_month(n) <= last,
                   EVIDownload.latest_per_month(names))
    expected = length(first:Month(1):last)
    length(names) == expected ||
        @warn "averaging $(length(names)) months where $expected were expected; run evi_download.jl first"
    println("Averaging ", length(names), " monthly files ...")
    flush(stdout)
    result = mean_evi(joinpath.(dir, names))
    write_mean(out, result.evi, result.months)
    write_provenance(out, source = "MOD13C2 v061 monthly files under data/evi",
                     product = "MODIS/Terra Vegetation Indices Monthly L3 Global 0.05Deg CMG V061",
                     doi = EVIDownload.DOI, citation = EVIDownload.CITATION,
                     first_month = string(EVIDownload.granule_month(names[1])),
                     last_month = string(EVIDownload.granule_month(names[end])),
                     months = length(names),
                     files = names,
                     reliability_kept = "0 (good) and 1 (marginal)",
                     cells_without_a_month = count(==(0), result.months))
    println("wrote ", out)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    EVIMean.main(ARGS)
end
