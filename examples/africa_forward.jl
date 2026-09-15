# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The forward run: the settled state `africa_burnin.jl` wrote, run from 1940 through the ERA5 months
# in calendar order, recording normalised subcommunity alpha diversity in every cell each 1 January.
# One replicate per seed; `AFRICA_NEUTRAL=true` repeats the 1940 climate throughout instead, the
# neutral control.
#
#     AFRICA_SCALE=small julia --project examples/africa_forward.jl
#     AFRICA_SEED=2 mpiexecjl --project -n 8 julia --project examples/africa_forward.jl
#
# | variable          | default        | meaning                                                  |
# |-------------------|----------------|----------------------------------------------------------|
# | `AFRICA_SCALE`    | `large`        | `small`: the same species and grid as a small burn-in, 1 year |
# | `AFRICA_CELLSIZE` | chosen         | as for the burn-in, which it must match                  |
# | `AFRICA_SEED`     | `1`            | the replicate                                            |
# | `AFRICA_NEUTRAL`  | `false`        | `true` to repeat the 1940 climate                        |
# | `AFRICA_SAVEDIR`  | `./Africa_run` | where the settled state is read and the diversity written |
# | `AFRICA_ENVELOPES`| the CEH file   | as for the burn-in                                       |

using Africa
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
const DATA = joinpath(pkgdir(Africa), "data")
const ENVELOPES = get(ENV, "AFRICA_ENVELOPES",
                      expanduser("~/git/plant_trait/data/global_plant_climate_envelopes.csv"))
const SEED = parse(Int, get(ENV, "AFRICA_SEED", "1"))
const NEUTRAL = get(ENV, "AFRICA_NEUTRAL", "false") == "true"
const SAVEDIR = get(ENV, "AFRICA_SAVEDIR", joinpath(pwd(), "Africa_run"))
const YEARS = smallscale() ? 1 : 2025 - 1940 + 1

say(xs...) = ISROOT && (println(xs...); flush(stdout))

# The same species, in the same order, as the burn-in built: the settled state is indexed by them.
envelopes = CSV.read(ENVELOPES, DataFrame,
                     select = [:species, :tmin, :tp, :ssr, :swvl1])
widths = CSV.read(joinpath(DATA, "species_widths.csv"), DataFrame)
pool = species_pool(envelopes, widths).pool
if smallscale()
    records = CSV.read(joinpath(DATA, "records_africa.csv"), DataFrame,
                       select = [:species])
    inpool = Set(pool.species)
    filter!(r -> r.species in inpool, records)
    counts = combine(groupby(records, :species), nrow => :n)
    chosen = Set(first(sort(counts, :n, rev = true), 200).species)
    pool = filter(r -> r.species in chosen, pool)
end
species = africa_species(pool, seed = SEED,
                         provenance = envelopes_record(ENVELOPES))
positioning = positioning_layers()
chosen = choose_cellsize(nrow(pool), regime = positioning.regime,
                         supply = positioning.supply)
area = africa_area(chosen.cellsize, regime = positioning.regime,
                   supply = positioning.supply,
                   verbosity = ISROOT ? :normal : :silent)

timing = NEUTRAL ? HeldClimateYear(1940) : ForwardClimate()
eco = africa_ecosystem(species, area, timing, seed = SEED)
continue_from!(eco, load_abundances(joinpath(SAVEDIR, "start00.jld2")),
               provenance = gbif_record("records_africa.csv"))
@time recorder = record_diversity!(eco, YEARS)

if ISROOT
    name = (NEUTRAL ? "neutral" : "forward") * "_seed$(SEED)"
    jldsave(joinpath(SAVEDIR, name * ".jld2"), diversity = recorder.storage,
            qs = [0.0, 1.0, 2.0, Inf], years = collect(1940:(1940 + YEARS)),
            active = vec(parent(eco.habitat.active)))
    write_provenance(joinpath(SAVEDIR, name * ".provenance.toml"),
                     EcoSISTEM.provenance(recorder))
    say("diversity written to ", joinpath(SAVEDIR, name * ".jld2"))
end
