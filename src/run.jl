# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Building and running the Africa ecosystem: the backwards burn-in that seeds it from the occurrence
# records, the years of held climate that settle it, and the forward run that records its
# diversity.

"""
    positioning_layers()

Return the `regime` and `supply` a study area over Africa is positioned and costed with: WorldClim's
annual mean temperature and a uniform sunlight. Reading the ERA5 files to decide the grid costs
minutes where these cost seconds, and the grid is the same projected lattice whichever positions
it; the run's own layers are then built onto it.
"""
function positioning_layers()
    return (regime = SourceSpec(WorldClim{BioClim}, :bio1),
            supply = UniformSpec(1.0kJ / (m^2 * Unitful.d),
                                 axis = SolarRadiation))
end

"""
    africa_ecosystem(species, area, timing = ForwardClimate(); seed = 1,
                     epoch = Dates.Date(1940, 1, 1), dir = era5_dir())

Build the ecosystem of `species` on `area` with the climate laid out as `timing` says, stepped by
`MeanMonths()` from `epoch`, serial or across MPI ranks as the process was launched.

# Arguments

  - `species`: the species list, from [`africa_species`](@ref).
  - `area`: the study area, from [`africa_area`](@ref).
  - `timing`: an [`AbstractClimateTiming`](@ref).
  - `seed`: seeds the run's random streams; replicates differ by seed.
  - `epoch`: the date elapsed time zero falls on.
  - `dir`: the directory of the ERA5 monthly files.
"""
function africa_ecosystem(species, area::StudyArea,
                          timing::AbstractClimateTiming = ForwardClimate();
                          seed::Integer = 1,
                          epoch::Dates.Date = Dates.Date(1940, 1, 1),
                          dir::AbstractString = era5_dir())
    layers = africa_layers(timing, dir = dir)
    habitat = GridHabitat(regime = layers.regime, supply = layers.supply,
                          area = area)
    return build_ecosystem(species, habitat, seed = seed, epoch = epoch,
                           calendar = MeanMonths(), distributed = :auto)
end

"""
    gbif_record(path)

Return the [`InputRecord`](@ref) of the GBIF occurrence download the seeding records come from,
read from `path`, for a run's provenance.
"""
function gbif_record(path::AbstractString)
    return EcoSISTEM.InputRecord(role = :abundance,
                                 dataset = "GBIF occurrence download",
                                 path = basename(path),
                                 doi = "10.15468/dl.n4a2sa",
                                 licence = "CC BY-NC 4.0",
                                 citation = "GBIF.org (2026) GBIF Occurrence Download " *
                                            "https://doi.org/10.15468/dl.n4a2sa")
end

# Each species' demand per individual on the resource named `name`, as a vector.
_demands(eco, name::Symbol) = getproperty(eco.spplist.demand, name).resource

# Each cell's supply of the resource named `name`, as a vector in cell order.
function _supplies(eco, name::Symbol)
    return vec(parent(getproperty(eco.habitat.supply, name).matrix))
end

"""
    burnin_table(eco, records, grid; density = 1.0 / m^2, start = (year = 2025, month = 12))

Return what a backwards burn-in of `eco` adds, as a named tuple: `table`, the abundance table for an
`AddAbundanceTable`, and the numbers that went into it - the scale `c`, the records `used`, those
`outside` the fine grid and those in an `inactive` cell, and the individuals `started` with.

The table's rows at elapsed zero fill each cell to what its resources hold at the run's first
month, taking species' shares from how many of their records the cell has; every record is then
added again at the month it was made, by its seeding factor.

# Arguments

  - `eco`: the ecosystem, built with [`ReversedClimate`](@ref).
  - `records`: the seeding records of `eco`'s species, as [`seeding_rows`](@ref) takes them.
  - `grid`: the fine seeding grid, from [`read_seeding_grid`](@ref).
  - `density`: plants per unit of available ground, as the species were built with.
  - `start`: the month the backwards run begins at, the last month of the climate.
"""
function burnin_table(eco, records::DataFrame, grid::NamedTuple;
                      density::Unitful.Quantity = 1.0 / m^2,
                      start::NamedTuple = (year = 2025, month = 12))
    active = vec(parent(eco.habitat.active))
    space = _supplies(eco, :space)
    landarea = sum(space[active])
    seeded = seeding_rows(records, grid, eco.habitat, landarea,
                          density = density, start = start)
    names = eco.spplist.names
    index = Dict(n => i for (i, n) in enumerate(names))
    shape = size(parent(eco.habitat.active))
    linear = LinearIndices(shape)
    rows = seeded.rows
    cells = [linear[y, x] for (y, x) in zip(rows.y, rows.x)]
    counts = SparseArrays.sparse([index[s] for s in rows.species], cells, 1.0,
                                 length(names), prod(shape))
    resources = (:space, :light, :water)
    filled = fill_to_capacity(counts, map(r -> _supplies(eco, r), resources),
                              map(r -> _demands(eco, r), resources))
    positions = CartesianIndices(shape)[filled.cell]
    table = build_abundance_table(eco,
                                  (species = vcat(names[filled.species],
                                                  rows.species),
                                   y = vcat([p[1] for p in positions], rows.y),
                                   x = vcat([p[2] for p in positions], rows.x),
                                   time = vcat(zeros(length(filled.cell)) .*
                                               1.0month_mean_duration,
                                               rows.time),
                                   count = vcat(Float64.(filled.count),
                                                rows.count)))
    return (table = table, c = seeded.c, used = seeded.used,
            outside = seeded.outside, inactive = seeded.inactive,
            started = sum(filled.count, init = 0))
end

"""
    run_burnin!(eco, table, months; provenance = EcoSISTEM.InputRecord[])

Run `eco`, built with [`ReversedClimate`](@ref), back through `months` months of climate in steps of
one `month_mean_duration`, adding the rows of `table` from [`burnin_table`](@ref) as their times
pass. `provenance` records the data the table was made from.
"""
function run_burnin!(eco, table, months::Integer;
                     provenance = EcoSISTEM.InputRecord[])
    seeding = Intervention(EveryStep(), ActiveCells(),
                           AddAbundanceTable(table), provenance = provenance)
    simulate!(eco, months * 1.0month_mean_duration, 1.0month_mean_duration,
              intervention = seeding)
    return eco
end

"""
    load_abundances(path)

Return the `species x cells` abundance matrix a `SaveAbundance` recorder wrote to `path`.
"""
load_abundances(path::AbstractString) = JLD2.load(path, "abun")

"""
    abundance_table(eco, abundances)

Return an abundance table placing every individual of `abundances`, a `species x cells` matrix
such as one saved from another ecosystem on the same grid and species, into `eco` at elapsed zero.
"""
function abundance_table(eco, abundances::AbstractMatrix{<:Integer})
    species, cells, counts = Int[], Int[], Int[]
    for (cell, column) in enumerate(eachcol(abundances))
        for (sp, n) in enumerate(column)
            n > 0 || continue
            push!(species, sp)
            push!(cells, cell)
            push!(counts, n)
        end
    end
    return (species = species, cell = cells, count = counts,
            time = zeros(length(cells)) .* 1.0month_mean_duration)
end

"""
    continue_from!(eco, abundances; provenance = EcoSISTEM.InputRecord[])

Set `eco`, freshly built, to hold `abundances`, a `species x cells` matrix saved from an ecosystem
on the same grid and species, by adding them all at elapsed zero. `provenance` records where the
abundances came from.
"""
function continue_from!(eco, abundances::AbstractMatrix{<:Integer};
                        provenance = EcoSISTEM.InputRecord[])
    EcoSISTEM.applyinterventions!(eco,
                                  Intervention(EveryStep(), ActiveCells(),
                                               AddAbundanceTable(abundance_table(eco,
                                                                                 abundances)),
                                               provenance = provenance),
                                  0.0month_mean_duration,
                                  1.0month_mean_duration,
                                  0)
    return eco
end

"""
    record_diversity!(eco, years; qs = [0.0, 1.0, 2.0, Inf])

Run `eco` forward `years` years in monthly steps and return the recorder holding its normalised
subcommunity alpha diversity in every cell at each order of `qs`, on every 1 January from the start:
`storage[:, :, k]` holds cells by orders at the `k`th.
"""
function record_diversity!(eco, years::Integer;
                           qs::AbstractVector{<:Real} = [0.0, 1.0, 2.0, Inf])
    cells = length(parent(eco.habitat.active))
    recorder = RecordDiversity(zeros(cells, length(qs), years + 1),
                               norm_sub_alpha, qs)
    simulate!(recorder, eco, years * 1.0year, 1.0month_mean_duration,
              every = EveryYear())
    return recorder
end
