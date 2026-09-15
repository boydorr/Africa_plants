# SPDX-License-Identifier: LGPL-3.0-or-later
#
# How many individuals each occurrence record seeds. EVI stands in for how many plants a place
# holds and a record count for how hard it was searched, so a record seeds more where the ground is
# greener or less visited. The scale is fixed so that the records of the whole continent would seed
# `density` plants on every square metre of available ground if every fine cell had been visited.

"""
    MIN_EVI_MONTHS

The fewest good or marginal months a 0.05 degree cell's mean EVI must be taken over to be trusted;
a cell with fewer has no EVI, and its records seed only themselves.
"""
const MIN_EVI_MONTHS = 36

"""
    seeding_scale(evi, fraction, cellarea, landarea; density = 1.0 / m^2)

Return `c`, the individuals a fine cell seeds per unit of `EVI x fraction x cellarea`, an inverse
area: `density` times `landarea`, divided by the sum of that product over the fine cells.

# Arguments

  - `evi`: each fine cell's mean EVI, `NaN` where it has none. A negative mean counts as no
    vegetation.
  - `fraction`: the share of each fine cell available to plants, 0 to 1.
  - `cellarea`: each fine cell's true area.
  - `landarea`: the available ground the seeds are meant to cover, the sum over the model grid's
    active cells of their available fraction times their area.
  - `density`: plants per unit area once seeded.

The three arrays share one shape and hold only the fine cells inside the simulated land; any cell
outside it must be `NaN` in `evi`.
"""
function seeding_scale(evi::AbstractArray{<:Real},
                       fraction::AbstractArray{<:Real},
                       cellarea::AbstractArray{<:Unitful.Area},
                       landarea::Unitful.Area;
                       density::Unitful.Quantity = 1.0 / m^2)
    size(evi) == size(fraction) == size(cellarea) ||
        error("`evi`, `fraction` and `cellarea` must share one shape; got $(size(evi)), " *
              "$(size(fraction)) and $(size(cellarea))")
    total = zero(float(eltype(cellarea)))
    for (e, f, a) in zip(evi, fraction, cellarea)
        isnan(e) || (total += max(e, 0) * f * a)
    end
    total > zero(total) ||
        error("no fine cell has vegetation to scale the seeding by")
    return uconvert(m^-2, density * landarea / total)
end

"""
    seeding_factors(c, evi, fraction, cellarea, records)

Return the individuals each record in a fine cell seeds, `max(c x EVI x fraction x cellarea /
records, 1)`, as an array of the cells' shape: a record always seeds at least the plant it saw, and
a cell without an EVI seeds only its records. Cells with no records hold `NaN`.

# Arguments

  - `c`: the scale from [`seeding_scale`](@ref).
  - `evi`, `fraction`, `cellarea`: as for [`seeding_scale`](@ref).
  - `records`: how many records each fine cell holds over the whole seeding period.
"""
function seeding_factors(c::Unitful.Quantity, evi::AbstractArray{<:Real},
                         fraction::AbstractArray{<:Real},
                         cellarea::AbstractArray{<:Unitful.Area},
                         records::AbstractArray{<:Integer})
    size(evi) == size(fraction) == size(cellarea) == size(records) ||
        error("`evi`, `fraction`, `cellarea` and `records` must share one shape")
    factors = fill(NaN, size(evi))
    for (i, (e, f, a, n)) in enumerate(zip(evi, fraction, cellarea, records))
        n > 0 || continue
        factors[i] = isnan(e) ? 1.0 :
                     max(uconvert(NoUnits, c * max(e, 0) * f * a) / n, 1.0)
    end
    return factors
end

"""
    months_back(year, month; start = (year = 2025, month = 12))

Return how many months each record's `(year, month)` lies before `start`, the month a backwards
run begins at: zero for a record in that month, negative for one after it.
"""
function months_back(year::AbstractVector{<:Integer},
                     month::AbstractVector{<:Integer};
                     start::NamedTuple = (year = 2025, month = 12))
    first = start.year * 12 + start.month
    return [first - (y * 12 + mo) for (y, mo) in zip(year, month)]
end

"""
    fine_cells(lat, lon, grid)

Return the linear index into the fine grid `grid` of the cell holding each position at latitudes
`lat` and longitudes `lon`, `0` for a position outside it. `grid` has the cell-centre vectors
`lon`, ascending, and `lat`, descending, and its arrays are indexed `[lon, lat]`.
"""
function fine_cells(lat::AbstractVector{<:Real}, lon::AbstractVector{<:Real},
                    grid::NamedTuple)
    step = grid.lon[2] - grid.lon[1]
    west = grid.lon[1] - step / 2
    north = grid.lat[1] + step / 2
    nx, ny = length(grid.lon), length(grid.lat)
    linear = LinearIndices((nx, ny))
    cells = zeros(Int, length(lat))
    for (k, (y, x)) in enumerate(zip(lat, lon))
        i = floor(Int, (x - west) / step) + 1
        j = floor(Int, (north - y) / step) + 1
        (1 <= i <= nx && 1 <= j <= ny) && (cells[k] = linear[i, j])
    end
    return cells
end

"""
    read_seeding_grid(path)

Return the fine seeding grid written by `data/src/seeding_grid.jl` as a named tuple of the cell
centres `lon` (ascending) and `lat` (descending) in degrees, and `evi`, `fraction` and `area`
indexed `[lon, lat]`: the mean EVI and the available share, `NaN` where missing, and each cell's
area.
"""
function read_seeding_grid(path::AbstractString)
    return NCDatasets.NCDataset(path) do ds
        lon = Float64.(ds["lon"][:])
        lat = Float64.(ds["lat"][:])
        evi = Float64.(coalesce.(ds["evi"][:, :], NaN))
        fraction = Float64.(coalesce.(ds["fraction"][:, :], NaN))
        area = repeat(reshape(Float64.(ds["area"][:]), 1, length(lat)),
                      length(lon), 1) .* m^2
        return (lon = lon, lat = lat, evi = evi, fraction = fraction,
                area = area)
    end
end

"""
    seeding_rows(records, grid, habitat, landarea; density = 1.0 / m^2,
                 start = (year = 2025, month = 12))

Return what a backwards run seeds from `records` onto `habitat`, as a named tuple: `rows`, the
columns `species`, `y`, `x`, `time` and `count` that `build_abundance_table` takes, one per record;
the scale `c` from [`seeding_scale`](@ref); and how many records were `used`, fell `outside` the
fine grid, or fell in no active cell of `habitat` (`inactive`).

A record's count is its fine cell's factor from [`seeding_factors`](@ref), with the scale normalised
over the fine cells whose centres lie in active cells of `habitat`, and its time is how many months
before `start` it was made, as an elapsed time in months of `month_mean_duration`.

# Arguments

  - `records`: one row per record, with `species`, `decimalLatitude`, `decimalLongitude`, `year`
    and `month`, holding only species `habitat`'s ecosystem simulates.
  - `grid`: the fine grid, as [`read_seeding_grid`](@ref) returns it.
  - `habitat`: the habitat the records are seeded into.
  - `landarea`: the available ground over `habitat`'s active cells.
  - `density`, `start`: plants per unit of available ground, and the month the run begins at.
"""
function seeding_rows(records::DataFrame, grid::NamedTuple,
                      habitat::GridHabitat,
                      landarea::Unitful.Area;
                      density::Unitful.Quantity = 1.0 / m^2,
                      start::NamedTuple = (year = 2025, month = 12))
    nx, ny = length(grid.lon), length(grid.lat)
    centrelat = vec(repeat(reshape(grid.lat, 1, ny), nx, 1))
    centrelon = vec(repeat(grid.lon, 1, ny))
    inland = reshape(grid_cells(habitat, centrelat, centrelon) .> 0,
                     nx, ny)
    fraction = replace(grid.fraction, NaN => 0.0)
    evi = copy(grid.evi)
    evi[.!inland] .= NaN
    c = seeding_scale(evi, fraction, grid.area, landarea, density = density)

    lat, lon = records.decimalLatitude, records.decimalLongitude
    fine = fine_cells(lat, lon, grid)
    counts = zeros(Int, nx, ny)
    for f in fine
        f > 0 && (counts[f] += 1)
    end
    factors = seeding_factors(c, grid.evi, fraction, grid.area, counts)
    model = grid_cells(habitat, lat, lon)
    keep = (fine .> 0) .& (model .> 0)
    positions = CartesianIndices(size(parent(habitat.active)))[model[keep]]
    back = months_back(records.year[keep], records.month[keep], start = start)
    rows = (species = records.species[keep],
            y = [p[1] for p in positions], x = [p[2] for p in positions],
            time = back .* 1.0month_mean_duration,
            count = factors[fine[keep]])
    return (rows = rows, c = c, used = count(keep),
            outside = count(==(0), fine),
            inactive = count((fine .> 0) .& (model .== 0)))
end
