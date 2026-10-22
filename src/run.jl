# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Building and running the Africa ecosystem: the backwards burn-in that seeds it from the occurrence
# records, the years of held climate that settle it, and the forward run that records its
# diversity.

"""
    positioning_layers()

Return the `regime` and `supply` a study area over Africa is positioned and costed with: no regime,
and the ground available to plants as the supply. Land cover has values over the sea as well as the
land, so it removes no coastal cell the outline touches, and the outline alone decides which cells
are simulated.
"""
function positioning_layers()
    return (regime = nothing, supply = available_ground())
end

"""
    africa_ecosystem(species, area, timing = ForwardClimate(); seed = 1,
                     epoch = Dates.Date(1940, 1, 1), dir = era5_dir(), distributed = :auto)

Build the ecosystem of `species` on `area` with the climate laid out as `timing` says, stepped by
`MeanMonths()` from `epoch`, serial or across MPI ranks as the process was launched.

# Arguments

  - `species`: the species list, from [`africa_species`](@ref).
  - `area`: the study area, from [`africa_area`](@ref).
  - `timing`: an [`AbstractClimateTiming`](@ref).
  - `seed`: seeds the run's random streams; replicates differ by seed.
  - `epoch`: the date elapsed time zero falls on.
  - `dir`: the directory of the ERA5 monthly files.
  - `distributed`: `:auto` builds a distributed ecosystem when the process is a live multi-rank MPI
    session and a serial one otherwise; `true` and `false` demand one or the other.
"""
function africa_ecosystem(species, area::StudyArea,
                          timing::AbstractClimateTiming = ForwardClimate();
                          seed::Integer = 1,
                          epoch::Dates.Date = Dates.Date(1940, 1, 1),
                          dir::AbstractString = era5_dir(), distributed = :auto)
    layers = africa_layers(timing, dir = dir)
    habitat = GridHabitat(regime = layers.regime, supply = layers.supply,
                          area = area)
    return build_ecosystem(species, habitat, seed = seed, epoch = epoch,
                           calendar = MeanMonths(), distributed = distributed)
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
    burnin_table(eco, records, grid; density = 1.0 / m^2, start = (year = 2025, month = 12),
                 halflife = 20.0year, snap = 3)

Return what a backwards burn-in of `eco` adds, as a named tuple: `table`, the abundance table for an
`AddAbundanceTable`, and the numbers that went into it - the EVI scale `c` and the continental
`scale`, the records `used`, those `outside` the fine grid, those `snapped` to a nearby active cell
and those left in an `inactive` cell, and the
individuals `started` with, `added` over the burn-in and added as `extra` from records after `start`.

Half the individuals are present at elapsed zero and half arrive as the run passes the months their
records were made, across the continent as a whole. The rows at elapsed zero fill each cell to what
its resources hold at the run's first month, the species' shares taken from their records in the
cell, each weighted by its seeding factor and by [`recency_weight`](@ref). Every record is then added
at its own month, by its seeding factor times one continental scale and at least one individual, the
scale chosen so that the records within the run add as many individuals as the start holds. Records
made after `start` are added at elapsed zero on top of both halves.

# Arguments

  - `eco`: the ecosystem, built with [`ReversedClimate`](@ref).
  - `records`: the seeding records of `eco`'s species, as [`seeding_rows`](@ref) takes them.
  - `grid`: the fine seeding grid, from [`read_seeding_grid`](@ref).
  - `density`: plants per unit of available ground, as the species were built with.
  - `start`: the month the backwards run begins at, the last month of the climate.
  - `halflife`: how long before `start` a record's weight in the starting shares halves.
  - `snap`: how many cells a record in an inactive cell may be moved to reach an active one.
"""
function burnin_table(eco, records::DataFrame, grid::NamedTuple;
                      density::Unitful.Quantity = 1.0 / m^2,
                      start::NamedTuple = (year = 2025, month = 12),
                      halflife::Unitful.Time = 20.0year, snap::Integer = 3)
    active = vec(parent(eco.habitat.active))
    space = _supplies(eco, :space)
    landarea = sum(space[active])
    seeded = seeding_rows(records, grid, eco.habitat, landarea,
                          density = density, start = start, snap = snap)
    names = eco.spplist.names
    index = Dict(n => i for (i, n) in enumerate(names))
    shape = size(parent(eco.habitat.active))
    linear = LinearIndices(shape)
    rows = seeded.rows
    cells = [linear[y, x] for (y, x) in zip(rows.y, rows.x)]
    shares = [max(f, 1.0) * recency_weight(t, halflife = halflife)
              for (f, t) in zip(rows.count, rows.time)]
    counts = SparseArrays.sparse([index[s] for s in rows.species], cells,
                                 shares, length(names), prod(shape))
    resources = (:space, :light, :water)
    filled = fill_to_capacity(counts, map(r -> _supplies(eco, r), resources),
                              map(r -> _demands(eco, r), resources))
    started = sum(filled.count, init = 0)
    within = [t >= zero(t) for t in rows.time]
    scaled = scale_to_total(rows.count[within], started)
    added = max.(scaled.scale .* rows.count, 1.0)
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
                                                added)))
    return (table = table, c = seeded.c, scale = scaled.scale,
            used = seeded.used, outside = seeded.outside,
            snapped = seeded.snapped, inactive = seeded.inactive,
            started = started,
            added = sum(added[within], init = 0.0),
            extra = sum(added[.!within], init = 0.0))
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
    DIVERSITY_ORDERS

The orders of diversity a run records: 0, the species richness, 1, 2 and infinity.
"""
const DIVERSITY_ORDERS = [0.0, 1.0, 2.0, Inf]

# Diversity.jl's metacommunity measures, by the name each is recorded under. Beta is the raw
# effective number of distinct subcommunities and rho its redundancy; the normalised pair weight
# each subcommunity by its share of the individuals.
const METACOMMUNITY_MEASURES = (raw_alpha = Diversity.raw_meta_alpha,
                                norm_alpha = Diversity.norm_meta_alpha,
                                raw_beta = Diversity.raw_meta_beta,
                                norm_beta = Diversity.norm_meta_beta,
                                raw_rho = Diversity.raw_meta_rho,
                                norm_rho = Diversity.norm_meta_rho,
                                gamma = Diversity.meta_gamma)

"""
    population(eco)

Return what `eco` holds as a named tuple of the `total` individuals, the `alive` species with any
individuals left, and the `largest` species' individuals. Under MPI each rank holds every cell of
its own species, so the counts are summed across the ranks and every rank must call it.
"""
function population(eco)
    perspecies = vec(sum(eco.abundances.matrix, dims = 2))
    return (total = round(Int, sum(perspecies)),
            alive = count(>(0), perspecies),
            largest = maximum(perspecies, init = 0))
end

function population(eco::MPIEcosystem)
    perspecies = vec(sum(eco.abundances.rows_matrix, dims = 2))
    comm = MPI.COMM_WORLD
    return (total = MPI.Allreduce(sum(perspecies), +, comm),
            alive = MPI.Allreduce(count(>(0), perspecies), +, comm),
            largest = MPI.Allreduce(maximum(perspecies, init = 0), max, comm))
end

"""
    abundance_totals(eco)

Return the individuals `eco` holds summed two ways, as a named tuple: `species`, each species' total
over the landscape, and `cells`, each cell's total over the species, in the orders of the species
list and of the grid's cells. They are small enough to keep for every recorded year, where the whole
landscape is not, and say which species a run lost and when, and where its individuals are. Under
MPI each rank holds every cell of its own species, so both are summed across the ranks and every
rank must call it.
"""
function abundance_totals(eco)
    matrix = eco.abundances.matrix
    return (species = vec(sum(matrix, dims = 2)),
            cells = vec(sum(matrix, dims = 1)))
end

function abundance_totals(eco::MPIEcosystem)
    rows = eco.abundances.rows_tuple
    matrix = eco.abundances.rows_matrix
    # Each rank's species go into their places in a landscape-wide vector, so that summing the
    # vectors across the ranks gives every rank the whole of it.
    species = zeros(Int, rows.total)
    species[rows.first:rows.last] = vec(sum(matrix, dims = 2))
    comm = MPI.COMM_WORLD
    return (species = MPI.Allreduce(species, +, comm),
            cells = MPI.Allreduce(vec(sum(matrix, dims = 1)), +, comm))
end

"""
    metacommunity_diversity(eco; qs = DIVERSITY_ORDERS)

Return the diversity of `eco` as a whole at each order in `qs`, as a named tuple of `raw_alpha`,
`norm_alpha`, `raw_beta`, `norm_beta`, `raw_rho`, `norm_rho` and `gamma`, each a vector over `qs`.
Every measure is Diversity.jl's own, over the landscape partitioned into its cells.

Under MPI no rank holds the whole landscape, so it is gathered first: every rank allocates a
`species x cells` array (3.8 GiB for the continent at 50 km), the root fills and measures it, and
the answers go back to every rank. Every rank must call this, and a run that records often is
better off in one process.
"""
function metacommunity_diversity(eco;
                                 qs::AbstractVector{<:Real} = DIVERSITY_ORDERS)
    return _measure(eco, qs)
end

# An ecosystem of its own is what Diversity.jl measures directly.
_measure(eco, qs) = map(m -> m(eco, qs).diversity, METACOMMUNITY_MEASURES)

# A distributed one has to be gathered first, and `gatherabundance` fills the array on the root
# alone, so only the root can measure it; the answers then go back to every rank, which keeps a
# caller that reports or records them the same on all of them.
function _measure(eco::MPIEcosystem, qs)
    gathered = EcoSISTEM.gatherabundance(eco)
    measured = MPI.Comm_rank(MPI.COMM_WORLD) == 0 ?
               _measure(Metacommunity(Float64.(gathered)), qs) :
               map(_ -> Float64[], METACOMMUNITY_MEASURES)
    return MPI.bcast(measured, 0, MPI.COMM_WORLD)
end

"""
    record_diversity(path, eco, stage, years; qs = DIVERSITY_ORDERS)

Append what `eco` holds to the table at `path`, one row per measure and order: the `stage` of the
run and how many `years` into it this is, the individuals, living species and largest species of
[`population`](@ref), and every measure of [`metacommunity_diversity`](@ref). Write the header first
if the file is new. Return the population and the measures, so that a caller can report them.

Under MPI every rank must call this, because the measures are gathered; only the root writes.
"""
function record_diversity(path::AbstractString, eco, stage::AbstractString,
                          years::Real;
                          qs::AbstractVector{<:Real} = DIVERSITY_ORDERS)
    counts = population(eco)
    measures = metacommunity_diversity(eco, qs = qs)
    if !MPI.Initialized() || MPI.Comm_rank(MPI.COMM_WORLD) == 0
        isfile(path) ||
            open(io -> println(io,
                               "stage,years,individuals,species,largest,measure,q,diversity"),
                 path, "w")
        open(path, "a") do io
            for (name, values) in pairs(measures), (q, value) in zip(qs,
                                                                     values)
                println(io, stage, ",", years, ",", counts.total, ",",
                        counts.alive, ",", counts.largest, ",", name, ",", q,
                        ",", value)
            end
        end
    end
    return (counts = counts, measures = measures)
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
