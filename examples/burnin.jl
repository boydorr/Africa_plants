using Africa
using Diversity
using JLD
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using JuliaDB
using JuliaDBMeta
using EcoSISTEM.ClimatePref
using Statistics
using EcoSISTEM
using Distributions
using AxisArrays
using Random

# Run burnin backwards from present day to 1901
function runSim()
    gbif = JuliaDB.load("data/GBIF_africa_fil")
    gbif = filter(g -> !ismissing(g.date), gbif)
    traits = JuliaDB.load("data/Africa_traits_fil")
    file = "data/Africa.tif"
    africa = readfile(file, -25°, 50°, -35°, 40°)[:, end:-1:1]

    # Set up grid
    numSpecies = length(traits); grd = (100,100); area = 64e6km^2;
    individuals = 0

    # Set up species requirements
    sizes = abs.(rand(Normal(0.01, 0.01), numSpecies))
    solarreq = collect(select(traits, :ssr))
    req1 = SolarRequirement(solarreq .* sizes  .* m^2)

    waterreq = collect(select(traits, :swvl1))
    req2 = VolWaterRequirement(waterreq .* sizes)

    req = ReqCollection2(req1, req2)

    tmean = collect(select(traits, :tmin_mean))
    tsd = collect(select(traits, :tmin_sd))
    tsd .+= (0.1 * maximum(tsd))
    temp_traits = GaussTrait(tmean, tsd)

    pmean = uconvert.(mm, collect(select(traits, :tp_mean)))
    psd = uconvert.(mm, collect(select(traits, :tp_sd)))
    psd .+= (0.1 * maximum(psd))
    prec_traits = GaussTrait(pmean, psd)

    av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)
    movement = BirthOnlyMovement(kernel, NoBoundary())

    abun = rand(Multinomial(individuals, numSpecies))
    trts = TraitCollection2(temp_traits, prec_traits)

    birth_rates = death_rates = fill(0.15, numSpecies) ./year
    param = PopGrowth{typeof(unit(birth_rates[1]))}(birth_rates, death_rates, 1.0, 1e-3, 1.0)

    native = fill(true, numSpecies)

    sppl = SpeciesList(numSpecies, trts, abun, req, movement, param, native)

    # Load CERA-20C/ERA climates
    JLD2.@load "data/Africa_temp.jld2"
    JLD2.@load "data/Africa_prec.jld2"
    JLD2.@load "data/Africa_water.jld2"
    JLD2.@load "data/Africa_solar.jld2"

    # Put together abiotic environment
    bud = BudgetCollection2(africa_solar, africa_water)
    active =  Array{Bool, 2}(.!isnan.(africa))
    ae1 = eraAE(africa_temp, africa_solar, active)
    ae2 = eraAE(africa_prec, africa_solar, active)
    hab = HabitatCollection2(ae1.habitat, ae2.habitat)
    ae = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, ae1.active, bud, ae1.names)

    rel1 = Gauss{eltype(ae.habitat.h1)}()
    rel2 = Gauss{eltype(ae.habitat.h2)}()
    rel = multiplicativeTR2(rel1, rel2)
    eco = Ecosystem(sppl, ae, rel)

    file = "data/Africa_effort.tif"
    effort = Array(readfile(file, -25°, 50°, -35°, 40°))[:, end:-1:1]
    effort[isnan.(effort)] .= 1
    effort[.!active] .= 0.0

    start = startingArray(gbif, numSpecies, true) .* hcat(effort[1:end]...)

    eco.abundances.matrix .+= round.(Int64, start)
    eco.abenv.budget.b2.matrix .+= 1e-9m^3
    fillarray = Array{Int64, 2}(undef, size(eco.abundances.matrix, 1), size(eco.abundances.matrix, 2))
    reverseBurnin!(eco, 118years - 1month, 1month, 10years, "../../sdb/Chapter5/burnin/", gbif, effort, fillarray)
    runburnin!(eco, 10years, 1month)
    return eco
end

output = runSim()

# Save output
diver = EcoSISTEM.SavedLandscape(output.abundances)
JLD2.@save "~/sdc/1901.jld2" diver

# Check number of species left
sum(mapslices(sum, output.abundances.matrix, dims =2)[:,1].>0)
