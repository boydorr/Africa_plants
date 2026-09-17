# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The layers a run is built from: the ERA5 monthly climate species are matched to, and the three
# resources they compete for, each confined to the ground in a cell available to plants.

"""
    AVAILABLE_CLASSES

The EarthEnv land-cover classes counted as ground plants can grow on: every class but open water,
snow and ice, and urban or built-up land.
"""
const AVAILABLE_CLASSES = (:needleleaf_trees, :evergreen_broadleaf_trees,
                           :deciduous_broadleaf_trees, :other_trees, :shrubs,
                           :herbaceous, :cultivated_and_managed,
                           :regularly_flooded, :barren)

"""
    era5_dir()

Return the directory the ERA5 monthly files are read from: `era5` under `AFRICA_DATADIR` if set,
and otherwise under the repository's `data/`.
"""
function era5_dir()
    return joinpath(get(ENV, "AFRICA_DATADIR",
                        joinpath(pkgdir(@__MODULE__), "data")), "era5")
end

"""
    available_ground()

Return the layer of ground available to plants in each cell, the summed shares of the
[`AVAILABLE_CLASSES`](@ref), as a space supply.
"""
function available_ground()
    return ConstructedRasterSpec(EarthEnv{LandCover},
                                 collect(AVAILABLE_CLASSES),
                                 axis = SurfaceArea) do bands...
        return sum(bands)
    end
end

"""
    AbstractClimateTiming

How a run's climate is laid out in time: [`ForwardClimate`](@ref), [`ReversedClimate`](@ref) or
[`HeldClimateYear`](@ref).
"""
abstract type AbstractClimateTiming end

"""
    ForwardClimate()

The ERA5 months in calendar order as a dated series, held at the last month once it is passed.
"""
struct ForwardClimate <: AbstractClimateTiming end

"""
    ReversedClimate()

The ERA5 months from the last back to the first, one `month_mean_duration` apart from elapsed
zero, as an undated series held at the first month once it is passed.
"""
struct ReversedClimate <: AbstractClimateTiming end

"""
    HeldClimateYear(year)

The twelve ERA5 months of `year` as a climatology of the months of the year, repeated for as long
as a run lasts.
"""
struct HeldClimateYear <: AbstractClimateTiming
    year::Int
end

# What a series of each timing does past its last slice, and what its coordinates mean.
_atend(::ForwardClimate) = HoldAtEnd()
_atend(::ReversedClimate) = HoldAtEnd()
_atend(::HeldClimateYear) = RepeatAtEnd()
_calendar(::ForwardClimate) = nothing
_calendar(::ReversedClimate) = UndatedSeries()
_calendar(::HeldClimateYear) = MonthOfYearSeries()

# The raster `r`, a monthly stack on real dates, laid out as `timing` says.
_retime(r, ::ForwardClimate) = r

function _retime(r, ::ReversedClimate)
    n = size(r.array, DimensionalData.Ti)
    return _slices(r, n:-1:1, (0:(n - 1)) .* 1.0month_mean_duration)
end

function _retime(r, timing::HeldClimateYear)
    dates = DimensionalData.lookup(r.array, DimensionalData.Ti)
    months = findall(t -> Dates.year(t) == timing.year, dates)
    length(months) == 12 ||
        error("the climate holds $(length(months)) months of $(timing.year), not twelve")
    return _slices(r, months, (1:12) .* 1.0month_mean_duration)
end

# The slices `which` of raster `r`, in that order, relabelled with elapsed `times`.
function _slices(r, which, times)
    array = r.array[DimensionalData.Ti(which)]
    return EcoSISTEM.ClimateRaster(ERA,
                                   DimensionalData.set(array,
                                                       DimensionalData.Ti =>
                                                           collect(times)),
                                   r.code)
end

# An ERA5 monthly variable read from `dir` as a dated series.
function _era5source(code::AbstractString, dir::AbstractString)
    return SourceSpec(ERA, code, directory = dir, atend = HoldAtEnd())
end

# An ERA5 monthly variable read from `dir`, laid out in time as `timing` says. Every timing goes
# through a combine, forward included, so that the three timings build onto the same active cells.
function _era5(code::AbstractString, dir::AbstractString,
               timing::AbstractClimateTiming)
    return ConstructedRasterSpec(_era5source(code, dir),
                                 axis = EcoSISTEM.layeraxis(ERA, code),
                                 atend = _atend(timing),
                                 calendar = _calendar(timing)) do values
        return _retime(values, timing)
    end
end

# A supply of ERA5 variable `code` on `axis`, confined to the ground available to plants and laid
# out as `timing` says. Values below zero, which regridding and ERA5's packing leave in dry cells,
# are no supply at all.
function _gated(code::AbstractString, axis::Type{<:EcoSISTEM.NicheAxis},
                dir::AbstractString, timing::AbstractClimateTiming)
    return ConstructedRasterSpec(_era5source(code, dir),
                                 EarthEnv{LandCover},
                                 collect(AVAILABLE_CLASSES),
                                 axis = axis, atend = _atend(timing),
                                 calendar = _calendar(timing)
                                 ) do values, bands...
        gated = max.(values, zero(eltype(values.array))) .* sum(bands)
        return _retime(gated, timing)
    end
end

"""
    africa_layers(timing = ForwardClimate(); dir = era5_dir())

Return the layers of a run as a named tuple of `regime`, the ERA5 2 m temperature and total
precipitation named `temperature` and `rainfall`, and `supply`, the available ground named `space`
and the ERA5 net solar radiation named `light`, the latter scaled by the share of the cell available
to plants. Water reaches species through the precipitation tolerance rather than as a resource: the
envelopes' soil-water figures rest on a reanalysis whose own publisher describes them as qualitative
(see the plan). The climate is laid out in time as
`timing`, an [`AbstractClimateTiming`](@ref), says; a dated one steps one month per
`month_mean_duration` in a run built with a `MeanMonths()` calendar.
"""
function africa_layers(timing::AbstractClimateTiming = ForwardClimate();
                       dir::AbstractString = era5_dir())
    return (regime = (temperature = _era5("t2m", dir, timing),
                      rainfall = _era5("tp", dir, timing)),
            supply = (space = available_ground(),
                      light = _gated("ssr", SolarRadiation, dir, timing)))
end
