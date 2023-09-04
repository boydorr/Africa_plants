using Africa
using Diversity
using JLD2
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using Statistics
using StatsBase
using EcoSISTEM
using Distributions
using AxisArrays
using Random
using DataFrames

# Run full simulation forwards in time 1901-present
function buildEco(start::EcoSISTEM.SavedLandscape, repeatYear::Bool, cacheFolder::String, cacheFile::String)
    JLD2.@load("data/Africa_traits.jld2")

    file = "data/Africa.tif"
    africa = readfile(file, -25.0°, 50.0°, -35.0°, 40.0°)[:, end:-1:1]

    # Set up grid
    numSpecies = nrow(traits_dat); grd = (100,100); area = 64e6km^2;
    individuals = 0

    # Set up species requirements
    solarreq = traits_dat.ssr
    req1 = SolarRequirement(solarreq .* m^2)

    waterreq = traits_dat.swvl1
    req2 = VolWaterRequirement(waterreq)

    req = ReqCollection2(req1, req2)

    tmean = traits_dat.tmin_mean
    tsd = traits_dat.tmin_sd
    tsd .+= (0.1 * maximum(tsd))
    temp_traits = GaussTrait(tmean, tsd)

    pmean = uconvert.(mm, traits_dat.tp_mean)
    psd = uconvert.(mm, traits_dat.tp_sd)
    psd .+= (0.1 * maximum(psd))
    prec_traits = GaussTrait(pmean, psd)

    av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)
    movement = BirthOnlyMovement(kernel, NoBoundary())

    abun = rand(Multinomial(individuals, numSpecies))
    trts = TraitCollection2(temp_traits, prec_traits)

    death_rates = abs.(rand(Normal(0.15, 0.001), numSpecies)) ./year
    birth_rates = death_rates
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

    eco.abundances = GridLandscape(start, (numSpecies, grd[1], grd[2]))
    divfuns = [norm_sub_alpha, sub_gamma]#, raw_sub_alpha, norm_sub_beta, raw_sub_beta, norm_sub_rho, raw_sub_rho, sub_gamma]
    q = [0.0]#, 1.0, 2.0, Inf]
    simDict = Dict("times" => 118years - 1month, "burnin" => 0year, "interval" => 12months, "timestep" => 1month, "divfuns" => divfuns, "q" => q, "cacheInterval" => 10years, "fileName" => cacheFile)
    lensim = length(0month:simDict["interval"]:simDict["times"])
    diver = zeros(length(divfuns), size(eco.abundances.matrix, 2), lensim, length(q))

    #runburnin!(eco, simDict["burnin"], simDict["timestep"])
    if repeatYear
        Africa.keepYear!(eco)
    end
    println("Starting simulation ...")
    runsim!(diver, eco, simDict, cacheFolder)
    return diver
end

# Load burnin data to start from
JLD2.@load "/home/claireh/sdc/1901.jld2"
# Run on workstation
output = buildEco(diver, false, "run/", "Africa_test")
# Save output
JLD2.@save "Africa_test.jld2" output

