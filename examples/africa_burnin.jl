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
# | `AFRICA_SCALE`    | `large`               | `small`: 200 species, the coarsest grid, 6 months back and 1 year held |
# | `AFRICA_CELLSIZE` | chosen                | `50km` to override the memory-based choice       |
# | `AFRICA_DISTRIBUTED` | `auto`             | `serial` or `mpi` to demand one                  |
# | `AFRICA_SEED`     | `1`                   | seeds the species and the run                    |
# | `AFRICA_HELD`     | `10`                  | years the first year's climate is repeated       |
# | `AFRICA_SAVEDIR`  | `./Africa_run`        | where the settled state is written               |
# | `AFRICA_ENVELOPES`| the CEH envelope file | `global_plant_climate_envelopes.csv`             |

using Africa
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using CSV
using DataFrames
using MPI

MPI.Initialized() || MPI.Init()
const ISROOT = MPI.Comm_rank(MPI.COMM_WORLD) == 0
const DATA = joinpath(pkgdir(Africa), "data")
const ENVELOPES = get(ENV, "AFRICA_ENVELOPES",
                      expanduser("~/git/plant_trait/data/global_plant_climate_envelopes.csv"))
const SEED = parse(Int, get(ENV, "AFRICA_SEED", "1"))
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

say(xs...) = ISROOT && (println(xs...); flush(stdout))

# The species and their records. A small run keeps the most-recorded species.
envelopes = CSV.read(ENVELOPES, DataFrame,
                     select = [:species, :tmin, :tp, :ssr, :swvl1])
widths = CSV.read(joinpath(DATA, "species_widths.csv"), DataFrame)
built = species_pool(envelopes, widths)
records = CSV.read(joinpath(DATA, "records_africa.csv"), DataFrame,
                   select = [:species, :decimalLatitude, :decimalLongitude,
                       :year, :month])
inpool = Set(built.pool.species)
filter!(r -> r.species in inpool, records)
pool = built.pool
if smallscale()
    counts = combine(groupby(records, :species), nrow => :n)
    chosen = Set(first(sort(counts, :n, rev = true), 200).species)
    pool = filter(r -> r.species in chosen, pool)
    filter!(r -> r.species in chosen, records)
end
say(nrow(pool), " species (", length(built.unmatched), " with no envelope, ",
    length(built.nowidth), " with no width), ", nrow(records), " records")

species = africa_species(pool, seed = SEED,
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
started = time()
@time simulate!(eco, MONTHS * 1.0month_mean_duration, 1.0month_mean_duration,
                every = EveryInterval(120 * 1.0month_mean_duration),
                intervention = seeding_rows) do occurrence
    say("burn-in: ", 10 * (occurrence.count - 1), " years back, ",
        round(Int, time() - started), " s elapsed")
    return nothing
end
SaveAbundance(SAVEDIR, "burnin")(eco, (count = 1,))

# The first year's climate repeated, from the state the backwards run reached.
held = africa_ecosystem(species, area, HeldClimateYear(1940), seed = SEED,
                        distributed = DISTRIBUTED)
backwards = load_abundances(joinpath(SAVEDIR, "burnin00.jld2"))
continue_from!(held, backwards,
               provenance = gbif_record("records_africa.csv"))
started = time()
@time simulate!(held, HELD * 1.0year, 1.0month_mean_duration,
                every = EveryInterval(12 * 1.0month_mean_duration)
                ) do occurrence
    say("held climate: year ", occurrence.count - 1, " of ", HELD, ", ",
        round(Int, time() - started), " s elapsed")
    return nothing
end
SaveAbundance(SAVEDIR, "start")(held, (count = 1,))
say("settled state written to ", joinpath(SAVEDIR, "start00.jld2"))
