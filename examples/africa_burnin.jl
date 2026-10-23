# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The backwards burn-in: every Africa species with an envelope and widths, filled from the occurrence
# records to what each cell's resources hold at the last month of the climate, run back through the
# ERA5 months to the first while each record is added at the month it was made, then settled for
# some years under the first year's climate repeated. Writes the settled state for
# `africa_forward.jl`.
#
#     AFRICA_SCALE=small julia --project examples/africa_burnin.jl
#     julia -t 8 --project examples/africa_burnin.jl
#     mpiexecjl --project -n 8 julia --project examples/africa_burnin.jl
#
# | variable          | default               | meaning                                          |
# |-------------------|-----------------------|--------------------------------------------------|
# | `AFRICA_SCALE`    | `large`               | `small`: 200 species, the coarsest grid, 6 months back and 1 year held; `medium`: 3,000 species across every level of recording, the coarsest grid, the full span |
# | `AFRICA_CELLSIZE` | chosen                | `50km` to override the memory-based choice       |
# | `AFRICA_DISTRIBUTED` | `auto`             | `serial` or `mpi` to demand one                  |
# | `AFRICA_SEED`     | `1`                   | seeds the run's random streams                   |
# | `AFRICA_SPECIES_SEED` | `1`               | seeds the species' own draws; the same in every run of an experiment |
# | `AFRICA_HELD`     | `10`                  | years the first year's climate is repeated       |
# | `AFRICA_SAVEDIR`  | `./Africa_run`        | where the settled state is written               |
# | `AFRICA_ENVELOPES`| the CEH envelope file | `global_plant_climate_envelopes.csv`             |

using AfricaPlants
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using CSV
using DataFrames
using JLD2
using MPI

MPI.Initialized() || MPI.Init()
const ISROOT = MPI.Comm_rank(MPI.COMM_WORLD) == 0
const DATA = joinpath(pkgdir(AfricaPlants), "data")
const ENVELOPES = get(ENV, "AFRICA_ENVELOPES",
                      expanduser("~/git/plant_trait/data/global_plant_climate_envelopes.csv"))
const SEED = parse(Int, get(ENV, "AFRICA_SEED", "1"))
# The species are the same in every replicate, so that every replicate can start from one settled
# state: their own draws take this seed and the run's random streams take `SEED`.
const SPECIES_SEED = parse(Int, get(ENV, "AFRICA_SPECIES_SEED", "1"))
# `auto` builds a distributed ecosystem under a live multi-rank MPI session and a serial one
# otherwise; `mpi` and `serial` demand one or the other.
const DISTRIBUTED = let want = get(ENV, "AFRICA_DISTRIBUTED", "auto")
    want == "auto" ? :auto :
    want == "mpi" ? true :
    want == "serial" ? false :
    error("AFRICA_DISTRIBUTED=$want is not one of auto, mpi, serial")
end
const SAVEDIR = get(ENV, "AFRICA_SAVEDIR", joinpath(pwd(), "Africa_run"))
const MONTHS = smallscale() ? 6 : 12 * (2025 - 1940 + 1)
const HELD = smallscale() ? 1 : parse(Int, get(ENV, "AFRICA_HELD", "10"))

# Where each stage's yearly record of the whole landscape goes: individuals, species and every
# metacommunity measure at each order, one file per stage so a rerun of one does not disturb another.
record(stage) = joinpath(SAVEDIR, "diversity_$(stage)_seed$(SEED).csv")

say(xs...) = ISROOT && (println(xs...); flush(stdout))

# One progress line from a record: where the run has got to, what it holds, and its gamma diversity.
function report(where, recorded, elapsed)
    counts = recorded.counts
    say(where, ", ", round(Int, elapsed), " s elapsed, ", counts.total,
        " individuals, ", counts.alive, " species, largest ",
        round(100 * counts.largest / counts.total, digits = 3), "%")
    return say("  gamma: ",
               join(round.(recorded.measures.gamma, digits = 1), ", "),
               " at q = ", join(DIVERSITY_ORDERS, ", "))
end

# The species and their records. A small run keeps the most-recorded species.
envelopes = CSV.read(ENVELOPES, DataFrame,
                     select = [:species, :genus, :family, :order, :tmin,
                         :tp, :ssr, :swvl1])
widths = CSV.read(joinpath(DATA, "species_widths.csv"), DataFrame)
built = species_pool(envelopes, widths)
records = CSV.read(joinpath(DATA, "records_africa.csv"), DataFrame,
                   select = [:species, :decimalLatitude, :decimalLongitude,
                       :year, :month])
inpool = Set(built.pool.species)
filter!(r -> r.species in inpool, records)
pool = scaled_pool(built.pool, records.species)
if nrow(pool) < nrow(built.pool)
    chosen = Set(pool.species)
    filter!(r -> r.species in chosen, records)
end
say(nrow(pool), " species (", length(built.unmatched),
    " with no envelope left out, ", length(built.nowidth),
    " on their relatives' width alone), ", nrow(records), " records")

species = africa_species(pool, seed = SPECIES_SEED,
                         provenance = envelopes_record(ENVELOPES))
positioning = positioning_layers()
chosen = choose_cellsize(nrow(pool), regime = positioning.regime,
                         supply = positioning.supply)
chosen.fits ||
    error("$(nrow(pool)) species at $(chosen.cellsize) do not fit in the memory available")
area = africa_area(chosen.cellsize, regime = positioning.regime,
                   supply = positioning.supply,
                   verbosity = ISROOT ? :normal : :silent)
say("grid ", chosen.cellsize)

# Backwards through the climate, seeding as the records' months pass.
eco = africa_ecosystem(species, area, ReversedClimate(), seed = SEED,
                       distributed = DISTRIBUTED)
seeding = burnin_table(eco, records,
                       read_seeding_grid(joinpath(DATA, "seeding_grid.nc")))
say("seeding: c = ", seeding.c, ", continental scale ", seeding.scale, "; ",
    seeding.used, " records used, ", seeding.inactive, " in no active cell; ",
    seeding.started, " individuals to start, ", seeding.added,
    " added over the burn-in, ", seeding.extra, " from records after it")
# Seeded as `run_burnin!` does, but through `simulate!` directly so that a run of many hours reports
# each decade it passes; the callback only prints, which is safe on every MPI rank.
seeding_rows = Intervention(EveryStep(), ActiveCells(),
                            AddAbundanceTable(seeding.table),
                            provenance = gbif_record("records_africa.csv"))
# Each species' and each cell's individuals in every recorded year of a stage, kept beside its
# diversity: which species a stage lost and when, and where the individuals are.
function totals(eco, years::Integer)
    return (species = zeros(Int, nrow(pool), years + 1),
            cells = zeros(Int, length(parent(eco.habitat.active)), years + 1))
end
function keep!(kept, eco, column::Integer)
    now = abundance_totals(eco)
    kept.species[:, column] = now.species
    kept.cells[:, column] = now.cells
    return kept
end
function save_totals(stage::AbstractString, kept)
    ISROOT &&
        jldsave(joinpath(SAVEDIR, "totals_$(stage)_seed$(SEED).jld2"),
                species = pool.species, species_totals = kept.species,
                cell_totals = kept.cells,
                years = collect(0:(size(kept.species, 2) - 1)))
    return nothing
end

backwards_totals = totals(eco, MONTHS ÷ 12)
started = time()
@time simulate!(eco, MONTHS * 1.0month_mean_duration, 1.0month_mean_duration,
                every = EveryInterval(12 * 1.0month_mean_duration),
                intervention = seeding_rows) do occurrence
    back = occurrence.count - 1
    keep!(backwards_totals, eco, occurrence.count)
    now = record_diversity(record("burnin"), eco, "burnin", back)
    back % 10 == 0 &&
        report("burn-in: $(back) years back", now, time() - started)
    return nothing
end
SaveAbundance(SAVEDIR, "burnin")(eco, (count = 1,))
save_totals("burnin", backwards_totals)

# The first year's climate repeated, from the state the backwards run reached.
held = africa_ecosystem(species, area, HeldClimateYear(1940), seed = SEED,
                        distributed = DISTRIBUTED)
backwards = load_abundances(joinpath(SAVEDIR, "burnin00.jld2"))
continue_from!(held, backwards,
               provenance = gbif_record("records_africa.csv"))
held_totals = totals(held, HELD)
started = time()
@time simulate!(held, HELD * 1.0year, 1.0month_mean_duration,
                every = EveryInterval(12 * 1.0month_mean_duration)
                ) do occurrence
    year = occurrence.count - 1
    keep!(held_totals, held, occurrence.count)
    now = record_diversity(record("hold"), held, "hold", year)
    report("held climate: year $(year) of $(HELD)", now, time() - started)
    return nothing
end
SaveAbundance(SAVEDIR, "start")(held, (count = 1,))
save_totals("hold", held_totals)
say("settled state written to ", joinpath(SAVEDIR, "start00.jld2"))
