# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The two variables species are matched to, each declared beneath the axis EcoSISTEM ships, so that
# it inherits the axis's unit and bounds and is reached by anything written for it, and each with a
# density width of its own. A suitability is the tolerance's density times that width, so the width
# sets the level every species pays at its optimum and not which species or cell does better.

# Each width is `sqrt(2 pi)` times the geometric mean of the pool's tolerance widths on the axis,
# 1.81 K and 0.965 mm/day, which puts a species of typical breadth at a suitability of one at its
# optimum: a narrower species scores above one there and a broader one below. It is a literal because
# the declaration pastes it into a method; `test/runtests.jl` holds it to the pool it came from,
# where that data is to hand, so it is recomputed when the widths change.

"""
    AnnualMinimumTemperature

The least monthly 2 m temperature of a year, the variable the species' temperature tolerances were
fitted on and the `temperature` regime holds. A [`Temperature`](@ref).
"""
@nicheaxis(AnnualMinimumTemperature <: Temperature, densitywidth = 4.5K)

"""
    AnnualMeanPrecipitation

The mean precipitation of a year, the variable the species' precipitation tolerances were fitted on
and the `rainfall` regime holds. A [`Precipitation`](@ref).
"""
@nicheaxis(AnnualMeanPrecipitation <: Precipitation,
           densitywidth = 2.4mm / Unitful.d)
