# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Where a geographic position falls on the projected simulation grid: occurrence records and the
# fine EVI cells are given in latitude and longitude, and the grid is in the run's projected CRS.

"""
    grid_cells(habitat, lat, lon; active = true, snap = 0)

Return the linear index of the cell of `habitat`'s grid holding each position at latitudes `lat`
and longitudes `lon` (degrees, WGS 84), `0` for a position off the grid or, with `active`, in an
inactive cell. The index is column-major over `(y, x)`, a column of the abundance matrix.

With `active` and a positive `snap`, a position in an inactive cell is moved to the nearest active
cell no more than `snap` cells away in either direction, nearest by distance between cell centres and,
among equally near cells, the first in column-major order; `0` where there is none.
"""
function grid_cells(habitat::GridHabitat, lat::AbstractVector{<:Real},
                    lon::AbstractVector{<:Real}; active::Bool = true,
                    snap::Integer = 0)
    length(lat) == length(lon) ||
        error("`lat` and `lon` must have one value per position")
    snap >= 0 || error("`snap` must not be negative")
    mask = parent(habitat.active)
    linear = LinearIndices(mask)
    places = [LatLong(y * °, x * °) for (y, x) in zip(lat, lon)]
    return map(EcoSISTEM.getcellat(habitat, places)) do cell
        isnothing(cell) && return 0
        (!active || mask[cell]) && return linear[cell]
        return _nearestactive(mask, cell, snap)
    end
end

# The linear index of the active cell of `mask` nearest `cell` within `radius` cells in each direction,
# by squared distance between centres, ties to the smaller linear index; `0` where none is active.
function _nearestactive(mask::AbstractMatrix{Bool}, cell::CartesianIndex{2},
                        radius::Integer)
    y, x = Tuple(cell)
    ny, nx = size(mask)
    linear = LinearIndices(mask)
    best, bestdistance = 0, typemax(Int)
    for k in 1:radius
        for yy in max(1, y - k):min(ny, y + k),
            xx in max(1, x - k):min(nx, x + k)
            max(abs(yy - y), abs(xx - x)) == k || continue
            mask[yy, xx] || continue
            distance = (yy - y)^2 + (xx - x)^2
            index = linear[yy, xx]
            if distance < bestdistance ||
               (distance == bestdistance && index < best)
                best, bestdistance = index, distance
            end
        end
        best > 0 && (k + 1)^2 > bestdistance && break
    end
    return best
end
