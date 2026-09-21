# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The forward run: the settled state `africa_burnin.jl` wrote, run from 1940 through the ERA5 months
# in calendar order, recording normalised subcommunity alpha diversity in every cell each 1 January,
# with each species' and each cell's individuals.
# One replicate per seed; `AFRICA_NEUTRAL=true` repeats the 1940 climate throughout instead, the
# neutral control. Every replicate starts from the one settled state and holds the same species,
# built from `AFRICA_SPECIES_SEED`; `AFRICA_SEED` seeds only the run's random streams. A control is
# the comparison for the forward run of the **same seed**, from which it differs only in climate. Each year both write the whole landscape's diversity to
# `diversity_<forward|neutral>_seed<N>.csv` beside the per-cell diversity they save at the end.
#
#     AFRICA_SCALE=small julia --project examples/africa_forward.jl
#     AFRICA_SEED=2 mpiexecjl --project -n 8 julia --project examples/africa_forward.jl
#
# | variable          | default        | meaning                                                  |
# |-------------------|----------------|----------------------------------------------------------|
# | `AFRICA_SCALE`    | `large`        | `small` or `medium`: the same species and grid as a burn-in at that scale; `small` runs 1 year |
# | `AFRICA_CELLSIZE` | chosen         | as for the burn-in, which it must match                  |
# | `AFRICA_DISTRIBUTED` | `auto`             | `serial` or `mpi` to demand one                  |
# | `AFRICA_SEED`     | `1`            | the replicate: seeds the run's random streams            |
# | `AFRICA_SPECIES_SEED` | `1`        | seeds the species' own draws, as the burn-in's did       |
# | `AFRICA_NEUTRAL`  | `false`        | `true` to repeat the 1940 climate                        |
# | `AFRICA_SAVEDIR`  | `./Africa_run` | where the settled state is read and the diversity written |
# | `AFRICA_START`    | `start00.jld2` | the settled state in `AFRICA_SAVEDIR` to start from       |
# | `AFRICA_SAVEEND`  | `true`         | `false` to skip saving the state the run ends in, which is as large as the settled state it started from |
# | `AFRICA_ENVELOPES`| the CEH file   | as for the burn-in                                       |

using Africa
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using CSV
using DataFrames
using JLD2
using DimensionalData: DimensionalData
using Diversity: norm_sub_alpha
using MPI

MPI.Initialized() || MPI.Init()
const ISROOT = MPI.Comm_rank(MPI.COMM_WORLD) == 0
const DATA = joinpath(pkgdir(Africa), "data")
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
const NEUTRAL = get(ENV, "AFRICA_NEUTRAL", "false") == "true"
const SAVEDIR = get(ENV, "AFRICA_SAVEDIR", joinpath(pwd(), "Africa_run"))
const START = get(ENV, "AFRICA_START", "start00.jld2")
const SAVEEND = get(ENV, "AFRICA_SAVEEND", "true") == "true"
const YEARS = smallscale() ? 1 : 2025 - 1940 + 1

say(xs...) = ISROOT && (println(xs...); flush(stdout))

# The same species, in the same order, as the burn-in built: the settled state is indexed by them.
envelopes = CSV.read(ENVELOPES, DataFrame,
                     select = [:species, :genus, :family, :order, :tmin,
                         :tp, :ssr, :swvl1])
widths = CSV.read(joinpath(DATA, "species_widths.csv"), DataFrame)
pool = species_pool(envelopes, widths).pool
if smallscale() || mediumscale()
    records = CSV.read(joinpath(DATA, "records_africa.csv"), DataFrame,
                       select = [:species])
    pool = scaled_pool(pool, records.species)
end
species = africa_species(pool, seed = SPECIES_SEED,
                         provenance = envelopes_record(ENVELOPES))
positioning = positioning_layers()
chosen = choose_cellsize(nrow(pool), regime = positioning.regime,
                         supply = positioning.supply)
area = africa_area(chosen.cellsize, regime = positioning.regime,
                   supply = positioning.supply,
                   verbosity = ISROOT ? :normal : :silent)

timing = NEUTRAL ? HeldClimateYear(1940) : ForwardClimate()
eco = africa_ecosystem(species, area, timing, seed = SEED,
                       distributed = DISTRIBUTED)
continue_from!(eco, load_abundances(joinpath(SAVEDIR, START)),
               provenance = gbif_record("records_africa.csv"))

# The recorder is called from a callback of our own rather than passed to `simulate!`, so that the
# whole landscape can be recorded and reported in the same place: EcoSISTEM records the year, we add
# the year's metacommunity measures to `diversity_<stage>_seed<N>.csv` and report each decade.
const QS = DIVERSITY_ORDERS
const STAGE = NEUTRAL ? "neutral" : "forward"
storage = zeros(length(parent(eco.habitat.active)), length(QS), YEARS + 1)
recorder = RecordDiversity(storage, norm_sub_alpha, QS)
# Each species' and each cell's individuals in every recorded year: which species were lost and
# when, and where the individuals are, without keeping the whole landscape.
species_totals = zeros(Int, nrow(pool), YEARS + 1)
cell_totals = zeros(Int, length(parent(eco.habitat.active)), YEARS + 1)
started = time()
@time simulate!(eco, YEARS * 1.0year, 1.0month_mean_duration,
                every = EveryYear()) do occurrence
    recorder(eco, occurrence)
    totals = abundance_totals(eco)
    species_totals[:, occurrence.count] = totals.species
    cell_totals[:, occurrence.count] = totals.cells
    year = 1940 + occurrence.count - 1
    now = record_diversity(joinpath(SAVEDIR,
                                    "diversity_$(STAGE)_seed$(SEED).csv"), eco,
                           STAGE, occurrence.count - 1)
    if occurrence.count == 1 || year % 10 == 0 || occurrence.count == YEARS + 1
        counts = now.counts
        say(year, ": ", round(Int, time() - started), " s elapsed, ",
            counts.total, " individuals, ", counts.alive, " species, largest ",
            round(100 * counts.largest / counts.total, digits = 3), "%")
        say("  gamma: ", join(round.(now.measures.gamma, digits = 1), ", "),
            " at q = ", join(QS, ", "))
    end
    return nothing
end

# The state a run ends in, so that another can carry on from it, as `AFRICA_START` names. It is as
# large as the settled state it started from, about 10 GB for the continent at 50 km, so a machine
# short of room turns it off.
if SAVEEND
    SaveAbundance(SAVEDIR, "$(STAGE)_seed$(SEED)_end")(eco, (count = 1,))
    say("final state written to ",
        joinpath(SAVEDIR, "$(STAGE)_seed$(SEED)_end00.jld2"))
end

if ISROOT
    name = (NEUTRAL ? "neutral" : "forward") * "_seed$(SEED)"
    # The grid's cell centres in metres, so the results can be drawn without rebuilding the area.
    function centres(dim)
        return ustrip.(m,
                       DimensionalData.Lookups.shiftlocus(DimensionalData.Lookups.Center(),
                                                          DimensionalData.lookup(eco.habitat.active,
                                                                                 dim)))
    end
    jldsave(joinpath(SAVEDIR, name * ".jld2"), diversity = recorder.storage,
            qs = QS, years = collect(1940:(1940 + YEARS)),
            species = pool.species, species_totals = species_totals,
            cell_totals = cell_totals,
            active = vec(parent(eco.habitat.active)),
            shape = size(parent(eco.habitat.active)),
            y = centres(DimensionalData.Y), x = centres(DimensionalData.X),
            crs = "EPSG:10592")
    write_provenance(joinpath(SAVEDIR, name * ".provenance.toml"),
                     EcoSISTEM.provenance(recorder))
    say("diversity written to ", joinpath(SAVEDIR, name * ".jld2"))
end
