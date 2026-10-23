# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The fine grid the seeding is scaled on: over the Africa box at the 0.05 degree resolution of the
# mean EVI layer, each cell's mean EVI (missing where fewer than `AfricaPlants.MIN_EVI_MONTHS` months went
# into it), the share of the cell available to plants, and its true area.
#
#     julia --project data/src/seeding_grid.jl [--force]
#
# The available share is the sum of EarthEnv's consensus land-cover classes other than open water,
# snow and ice and urban or built-up land, each a percentage of its 30 arcsecond cell, averaged over
# the six by six such cells making up each 0.05 degree cell; a cell with no land-cover value has no
# share. Reads `evi/evi_mean_2001_2025.nc` and the EarthEnv bands (downloaded by EcoSISTEM on first
# use); writes `seeding_grid.nc`.

module SeedingGrid

using AfricaPlants
using ArchGDAL
using EcoSISTEM
using NCDatasets
using RasterDataSources

include("common.jl")

const CELL = 0.05
const LANDCOVER_CELLS = 6
const LANDCOVER_NODATA = 0xff
const EARTH_RADIUS = 6_371_007.2
const AVAILABLE = (:needleleaf_trees, :evergreen_broadleaf_trees,
                   :deciduous_broadleaf_trees, :other_trees, :shrubs,
                   :herbaceous, :cultivated_and_managed, :regularly_flooded,
                   :barren)

"""
    box_indices(box; cell = CELL, west = -180.0, north = 90.0)

Return the column and row ranges of a global north-up grid of `cell` degrees, whose first cell's
corner is at `west` and `north`, covering the geographic `box` (with `north`, `west`, `south` and
`east` edges).
"""
function box_indices(box::NamedTuple; cell::Real = CELL, west::Real = -180.0,
                     north::Real = 90.0)
    columns = (round(Int, (box.west - west) / cell) + 1):round(Int,
                                                               (box.east - west) / cell)
    rows = (round(Int, (north - box.north) / cell) + 1):round(Int,
                                                              (north - box.south) / cell)
    return (columns = columns, rows = rows)
end

"""
    cell_areas(latitudes; cell = CELL)

Return the area in square metres of a `cell`-degree cell centred at each of `latitudes`, on the
sphere of the WGS 84 authalic radius.
"""
function cell_areas(latitudes::AbstractVector{<:Real}; cell::Real = CELL)
    half = cell / 2
    return [EARTH_RADIUS^2 * deg2rad(cell) *
            (sind(lat + half) - sind(lat - half)) for lat in latitudes]
end

"""
    block_fraction!(total, count, band, factor)

Add one land-cover band, percentages over cells `factor` times finer than `total`, to the running
per-cell sums: into `total` the band's values over each block divided by 100, and into `count` the
number of cells with a value, the last band's count being the one kept.
"""
function block_fraction!(total::AbstractMatrix{Float64},
                         count::AbstractMatrix{Int},
                         band::AbstractMatrix{UInt8},
                         factor::Integer)
    size(band) == size(total) .* factor ||
        error("the band is not $factor times finer than the grid")
    fill!(count, 0)
    @inbounds for j in axes(band, 2), i in axes(band, 1)
        v = band[i, j]
        v == LANDCOVER_NODATA && continue
        bi, bj = (i - 1) ÷ factor + 1, (j - 1) ÷ factor + 1
        total[bi, bj] += v / 100
        count[bi, bj] += 1
    end
    return total, count
end

# One EarthEnv band's window of the box, `width x height` with the first row at the north.
function _readband(path::AbstractString, box::NamedTuple)
    idx = box_indices(box, cell = CELL / LANDCOVER_CELLS)
    return ArchGDAL.read(path) do ds
        gt = ArchGDAL.getgeotransform(ds)
        isapprox(gt[1], -180.0) && isapprox(gt[4], 90.0) &&
        isapprox(gt[2], CELL / LANDCOVER_CELLS) ||
            error("$(basename(path)) is not the 30 arcsecond global grid: geotransform $gt")
        return ArchGDAL.read(ArchGDAL.getband(ds, 1), idx.rows, idx.columns)
    end
end

"""
    available_fraction(paths, box)

Return the share of each 0.05 degree cell of `box` available to plants, `width x height` with the
first row at the north: the land-cover bands at `paths` summed and averaged over the cells with a
value, `NaN` where none has one.
"""
function available_fraction(paths::AbstractVector{<:AbstractString},
                            box::NamedTuple)
    idx = box_indices(box)
    total = zeros(length(idx.columns), length(idx.rows))
    count = zeros(Int, size(total))
    for path in paths
        block_fraction!(total, count, _readband(path, box), LANDCOVER_CELLS)
    end
    return [n == 0 ? NaN : t / n for (t, n) in zip(total, count)]
end

function main(argv)
    out = joinpath(datadir(), "seeding_grid.nc")
    if isfile(out) && !forced(argv)
        println("Already have ", out,
                "; nothing to do (use --force to rebuild).")
        return nothing
    end
    box = AFRICA_BOX
    idx = box_indices(box)
    evipath = joinpath(datadir(), "evi", "evi_mean_2001_2025.nc")
    evi, months = NCDataset(evipath) do ds
        return (Float64.(coalesce.(ds["evi"][idx.columns, idx.rows], NaN)),
                Int.(ds["months"][idx.columns, idx.rows]))
    end
    evi[months .< AfricaPlants.MIN_EVI_MONTHS] .= NaN
    paths = [getraster(EarthEnv{LandCover}, class) for class in AVAILABLE]
    fraction = available_fraction(paths, box)
    lon = -180 .+ CELL .* (idx.columns .- 0.5)
    lat = 90 .- CELL .* (idx.rows .- 0.5)
    NCDataset(out, "c") do ds
        defDim(ds, "lon", length(lon))
        defDim(ds, "lat", length(lat))
        defVar(ds, "lon", lon, ("lon",), attrib = ["units" => "degrees_east"])
        defVar(ds, "lat", lat, ("lat",), attrib = ["units" => "degrees_north"])
        defVar(ds, "evi", Float32.(evi), ("lon", "lat"), fillvalue = NaN32,
               deflatelevel = 4,
               attrib = ["long_name" => "mean EVI, missing below $(AfricaPlants.MIN_EVI_MONTHS) months"])
        defVar(ds, "months", Int16.(months), ("lon", "lat"), deflatelevel = 4)
        defVar(ds, "fraction", Float32.(fraction), ("lon", "lat"),
               fillvalue = NaN32, deflatelevel = 4,
               attrib = ["long_name" => "share of the cell available to plants"])
        defVar(ds, "area", cell_areas(lat), ("lat",),
               attrib = ["units" => "m2"])
        return nothing
    end
    write_sidecar(out,
                  source = "evi/evi_mean_2001_2025.nc and EarthEnv consensus land cover",
                  landcover = "EarthEnv{LandCover}, classes " *
                              join(AVAILABLE, ", "),
                  landcover_doi = "10.1111/geb.12182",
                  min_evi_months = AfricaPlants.MIN_EVI_MONTHS,
                  cells = length(evi),
                  cells_with_evi = count(!isnan, evi),
                  cells_with_fraction = count(!isnan, fraction))
    println("wrote ", out)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    SeedingGrid.main(ARGS)
end
