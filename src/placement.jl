# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Where a geographic position falls on the projected simulation grid: occurrence records and the
# fine EVI cells are given in latitude and longitude, and the grid is in the run's projected CRS.

"""
    grid_cells(habitat, lat, lon; active = true)

Return the linear index of the cell of `habitat`'s grid holding each position at latitudes `lat`
and longitudes `lon` (degrees, WGS 84), `0` for a position off the grid or, with `active`, in an
inactive cell. The index is column-major over `(y, x)`, a column of the abundance matrix.
"""
function grid_cells(habitat::GridHabitat, lat::AbstractVector{<:Real},
                    lon::AbstractVector{<:Real}; active::Bool = true)
    length(lat) == length(lon) ||
        error("`lat` and `lon` must have one value per position")
    mask = parent(habitat.active)
    linear = LinearIndices(mask)
    places = [LatLong(y * °, x * °) for (y, x) in zip(lat, lon)]
    return map(EcoSISTEM.getcellat(habitat, places)) do cell
        return isnothing(cell) || (active && !mask[cell]) ? 0 : linear[cell]
    end
end
