# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The layers a run is built from: the ERA5 monthly climate species are matched to, and the three
# resources they compete for - ground, light and water - each confined to the ground in a cell
# available to plants.

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

"""
    AbstractRegimeStatistic

What a regime layer holds of a monthly climate variable: each [`MonthlyValue`](@ref), or the
[`TrailingMinimum`](@ref) or [`TrailingMean`](@ref) over the months ending in each.
"""
abstract type AbstractRegimeStatistic end

"""
    MonthlyValue()

Each month's own value.
"""
struct MonthlyValue <: AbstractRegimeStatistic end

"""
    TrailingMinimum(months = 12)

The least monthly value over the `months` ending in each month, the statistic the species'
temperature tolerances were fitted on.
"""
struct TrailingMinimum <: AbstractRegimeStatistic
    months::Int
    TrailingMinimum(months::Integer = 12) = new(months)
end

"""
    TrailingMean(months = 12)

The mean monthly value over the `months` ending in each month, the statistic the species'
precipitation tolerances were fitted on.
"""
struct TrailingMean <: AbstractRegimeStatistic
    months::Int
    TrailingMean(months::Integer = 12) = new(months)
end

# `statistic` over the months of `window`, which run along dimension `t`.
_reduced(::TrailingMinimum, window, t::Int) = minimum(window, dims = t)
_reduced(::TrailingMean, window, t::Int) = mean(window, dims = t)

# The raster `r`, a monthly stack in time order, with every month replaced by `statistic` over the
# months ending in it. A month with too few before it takes the first full run of months instead, so
# that none is summarised over a shorter stretch than the rest; a stack no longer than the run is
# therefore the same throughout.
_summarised(r, ::MonthlyValue) = r
function _summarised(r, statistic::AbstractRegimeStatistic)
    array = r.array
    t = DimensionalData.dimnum(array, DimensionalData.Ti)
    # Reduced as bare numbers and given their unit back: a minimum over quantities cannot start from
    # a slice holding `NaN`, which every slice with sea in it does.
    unit = Unitful.unit(eltype(array))
    data = ustrip.(unit, parent(array))
    span = min(statistic.months, size(data, t))
    out = similar(data)
    for k in axes(data, t)
        last = max(k, first(axes(data, t)) + span - 1)
        window = selectdim(data, t, (last - span + 1):last)
        selectdim(out, t, k) .= dropdims(_reduced(statistic, window, t),
                                         dims = t)
    end
    return EcoSISTEM.ClimateRaster(ERA,
                                   DimensionalData.rebuild(array,
                                                           data = out .* unit),
                                   r.code)
end

# The raster `r` summarised by `statistic` and laid out as `timing` says. A calendar order is
# summarised before it is reversed, so that a month looks back over the same months whichever way
# the run passes it; a held year is summarised over its own twelve months, which makes it constant.
function _laidout(r, statistic::AbstractRegimeStatistic,
                  timing::AbstractClimateTiming)
    return _retime(_summarised(r, statistic), timing)
end
function _laidout(r, statistic::AbstractRegimeStatistic,
                  timing::HeldClimateYear)
    return _summarised(_retime(r, timing), statistic)
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

# An ERA5 monthly variable read from `dir` onto `axis`, the axis the species' tolerance names,
# summarised by `statistic` and laid out in time as `timing` says. Every timing goes through a combine, forward included, so that the three timings
# build onto the same active cells.
function _era5(code::AbstractString, axis::Type{<:EcoSISTEM.NicheAxis},
               dir::AbstractString, timing::AbstractClimateTiming,
               statistic::AbstractRegimeStatistic)
    return ConstructedRasterSpec(_era5source(code, dir), axis = axis,
                                 atend = _atend(timing),
                                 calendar = _calendar(timing)) do values
        return _laidout(values, statistic, timing)
    end
end

# A supply of ERA5 variable `code` on `axis`, confined to the ground available to plants, summarised
# by `statistic` and laid out as `timing` says. The share of each cell that is available is the
# ground layer itself, [`available_ground`](@ref), so that every supply agrees with it about what
# ground there is. Values below zero, which regridding and ERA5's packing leave in dry cells, are no
# supply at all.
function _gated(code::AbstractString, axis::Type{<:EcoSISTEM.NicheAxis},
                dir::AbstractString, timing::AbstractClimateTiming,
                statistic::AbstractRegimeStatistic = MonthlyValue())
    return ConstructedRasterSpec(_era5source(code, dir), available_ground(),
                                 axis = axis, atend = _atend(timing),
                                 calendar = _calendar(timing)) do values, ground
        gated = max.(values, zero(eltype(values.array))) .* ground
        return _laidout(gated, statistic, timing)
    end
end

"""
    africa_layers(timing = ForwardClimate(); dir = era5_dir(),
                  temperature = TrailingMinimum(), rainfall = TrailingMean(),
                  water = TrailingMean())

Return the layers of a run as a named tuple of `regime`, the ERA5 2 m temperature and total
precipitation named `temperature` and `rainfall`, and `supply`: the available ground named `space`,
and the ERA5 net solar radiation named `light` and total precipitation named `water`, each scaled by
the share of the cell available to plants.

Each regime layer holds the statistic the species' tolerances were fitted on, so that a species is
matched against the same quantity its tolerance describes: the least monthly temperature and the
mean precipitation over the twelve months ending in each month. Rain is both: as a regime it says
which species a place suits, and as a supply how many plants the place can hold, so that a desert
cannot fill to the density of a forest. The water supply is the year's mean too, since a plant draws
on what the ground has stored, and a month with no rain at all would otherwise be a month with no
supply. Light stays monthly.

# Arguments

  - `timing`: an [`AbstractClimateTiming`](@ref), how the climate is laid out in time; a dated one
    steps one month per `month_mean_duration` in a run built with a `MeanMonths()` calendar.
  - `dir`: the directory of the ERA5 monthly files.
  - `temperature`: the [`AbstractRegimeStatistic`](@ref) the temperature layer holds.
  - `rainfall`: the statistic the precipitation layer holds.
  - `water`: the statistic the precipitation supply holds.
"""
function africa_layers(timing::AbstractClimateTiming = ForwardClimate();
                       dir::AbstractString = era5_dir(),
                       temperature::AbstractRegimeStatistic = TrailingMinimum(),
                       rainfall::AbstractRegimeStatistic = TrailingMean(),
                       water::AbstractRegimeStatistic = TrailingMean())
    return (regime = (temperature = _era5("t2m", AnnualMinimumTemperature, dir,
                                          timing, temperature),
                      rainfall = _era5("tp", AnnualMeanPrecipitation, dir,
                                       timing, rainfall)),
            supply = (space = available_ground(),
                      light = _gated("ssr", SolarRadiation, dir, timing),
                      water = _gated("tp", Precipitation, dir, timing, water)))
end
