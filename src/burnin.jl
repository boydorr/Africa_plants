# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The backwards burn-in: the continent filled from the occurrence records to what its resources
# hold, then run back through the climate while every record is added at the time it was made.

"""
    fill_to_capacity(counts, supply, demand)

Return the starting abundances as a named tuple of vectors `species`, `cell` and `count`, one per
species present in a cell: each cell's `counts` scaled so that the cell holds as many individuals as
its scarcest resource supports at the mean demand of the individuals seeded there, every species
keeping its share, rounded to whole individuals. A cell with no counts stays empty.

# Arguments

  - `counts`: a sparse `species x cells` matrix of individuals to take the shares from.
  - `supply`: one vector per resource of each cell's supply.
  - `demand`: one vector per resource, in the order of `supply`, of each species' demand per
    individual, in a unit dividing that resource's supply to a number.
"""
function fill_to_capacity(counts::SparseArrays.SparseMatrixCSC{<:Real},
                          supply::Tuple, demand::Tuple)
    length(supply) == length(demand) ||
        error("`supply` and `demand` must name the same resources")
    all(length(s) == size(counts, 2) for s in supply) ||
        error("each supply needs one value per cell")
    all(length(d) == size(counts, 1) for d in demand) ||
        error("each demand needs one value per species")
    rows = SparseArrays.rowvals(counts)
    values = SparseArrays.nonzeros(counts)
    species, cells, filled = Int[], Int[], Int[]
    for cell in axes(counts, 2)
        range = SparseArrays.nzrange(counts, cell)
        total = sum(view(values, range), init = 0.0)
        total > 0 || continue
        capacity = minimum(zip(supply, demand)) do (s, d)
            mean = sum(values[k] * d[rows[k]] for k in range) / total
            return uconvert(NoUnits, s[cell] / mean)
        end
        for k in range
            n = round(Int, values[k] * capacity / total)
            n > 0 || continue
            push!(species, rows[k])
            push!(cells, cell)
            push!(filled, n)
        end
    end
    return (species = species, cell = cells, count = filled)
end
