using Africa_plants
using Diversity
using JLD
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using JuliaDB
using JuliaDBMeta
using EcoSISTEM.ClimatePref
using Statistics
using StatsBase
using EcoSISTEM
using Distributions
using AxisArrays
using Random

function buildEco(outputfile::String, repeatYear::Bool, cacheFolder::String, cacheFile::String)
    gbif = JuliaDB.load("data/GBIF_africa_fil")
    gbif = filter(g -> !ismissing(g.date), gbif)
    traits = JuliaDB.load("Africa_traits_fil")

    file = "data/Africa.tif"
    africa = readfile(file, -25.0°, 50.0°, -35.0°, 40.0°)[:, end:-1:1]

    # Set up grid
    numSpecies = length(traits); grd = (100,100); area = 64e6km^2;
    individuals = 0

    # Set up species requirements
    solarreq = collect(select(traits, :ssr))
    req1 = SolarRequirement(solarreq .* m^2)

    waterreq = collect(select(traits, :swvl1))
    req2 = VolWaterRequirement(waterreq)

    req = ReqCollection2(req1, req2)

    tmean = collect(select(traits, :tmin_mean))
    tsd = collect(select(traits, :tmin_sd))
    tsd .+= (0.1 * maximum(tsd))
    temp_traits = GaussTrait(tmean, tsd)

    pmean = collect(select(traits, :tp_mean))
    psd = collect(select(traits, :tp_sd))
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

    ## Load CERA climates
    dir1 = "data/CERA"
    tempax1 = readCERA(dir1, "cera_20c_temp2m", "t2m")
    precax1 = readCERA(dir1, "cera_20c_totalprec", "tp")
    soilwaterax1 = readCERA(dir1, "cera_20c_soilwater1", "swvl1")
    solarradax1 = readCERA(dir1, "cera_20c_surfacenetsolar", "ssr")

    times = [collect(1979year:1month:(1980year - 1.0month)),
        collect(1980year:1month:(1990year - 1.0month)),
        collect(1990year:1month:(2000year - 1.0month)),
        collect(2000year:1month:(2010year - 1.0month)),
        collect(2010year:1month:(2019year- 1.0month))]

    # Load CERA-20C/ERA climates
    JLD2.@load "Africa_temp"
    JLD2.@load "Africa_prec"
    JLD2.@load "Africa_water"
    JLD2.@load "Africa_solar"

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

    output = JLD.load(outputfile, "diver")
    eco.abundances = GridLandscape(output, (numSpecies, grd[1], grd[2]))
    divfuns = [norm_sub_alpha, sub_gamma]#, raw_sub_alpha, norm_sub_beta, raw_sub_beta, norm_sub_rho, raw_sub_rho, sub_gamma]
    q = [0.0]#, 1.0, 2.0, Inf]
    simDict = Dict("times" => 118years - 1month, "burnin" => 0year, "interval" => 1year, "timestep" => 1month, "divfuns" => divfuns, "q" => q, "cacheInterval" => 10years, "fileName" => cacheFile)
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

# Run on workstation
for i in 1:10
    output = buildEco("Africa_diversity_1.jld", false, "run/", "Africa_$i")
    JLD.save("Africa_full_$i.jld", "diver", output)

    output = buildEco("Africa_diversity_1.jld", true, "neutral/", "Africa_$i")
    JLD.save("Africa_neutral_$i.jld", "diver", output)
    print("$i \n")
end
