include("../GIT/Chapter5/src/Chapter5.jl")
#include("../src/Chapter5.jl")
using .Chapter5
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
#using TimerOutputs
using Random
#using Plots
#pyplot()
import EcoSISTEM: _getdimension, _traitfun, AbstractAbiotic, AbstractTraitRelationship
function ecopopulate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AB, rel::R) where {AB <: AbstractAbiotic, R <: AbstractTraitRelationship}
  # Calculate size of habitat
  dim = _getdimension(abenv.habitat.h1)
  numsquares = dim[1] * dim[2]
  numspp = length(spplist.names)
  probabilities = [_traitfun(abenv.habitat.h1, spplist.traits.t1, rel.tr1, i, spp).* _traitfun(abenv.habitat.h2, spplist.traits.t2, rel.tr2, i, spp) for i in 1:numsquares, spp in 1:numspp]
  # Loop through species
  for i in eachindex(spplist.abun)
      if spplist.native[i]
        # Get abundance of species
        abun = rand(Multinomial(spplist.abun[i], probabilities[:, i]./sum(probabilities[:, i])))
        # Add individual to this location
        ml.matrix[i, :] .+= abun
     end
   end
end

function runSim()
    gbif = JuliaDB.load("data/Full_GBIF_africa")
    traits = JuliaDB.load("data/Africa_traits")
    traits = filter(tr -> (tr.swvl1 > 0.0m^3) & (tr.ssr > 0.0J/m^2) & (tr.prec_sd > 0.0mm) & (tr.tmin_sd > 0.0K), traits)
    gbif = filter(gb -> gb.SppID in select(traits, :SppID), gbif)

    file = "data/Africa.tif"
    africa = readfile(file, -25°, 50°, -35°, 40°)[:, end:-1:1]
    #heatmap(africa)

    # Set up grid
    numSpecies = length(traits); grd = (100,100); area = 64e6km^2;
    individuals = Int64(3.5e13)

    # Set up species requirements
    solarreq = collect(select(traits, :ssr))
    req1 = SolarRequirement(solarreq .* m^2)

    waterreq = collect(select(traits, :swvl1))
    req2 = VolWaterRequirement(waterreq)

    req = ReqCollection2(req1, req2)

    tmean = collect(select(traits, :tmin_mean))
    tsd = collect(select(traits, :tmin_sd))
    temp_traits = GaussTrait(tmean, tsd)

    pmean = collect(select(traits, :prec_mean))
    psd = collect(select(traits, :prec_sd))
    prec_traits = GaussTrait(pmean, psd)

    av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
    kernel = GaussianKernel(av_dist, 10e-10)
    movement = BirthOnlyMovement(kernel, NoBoundary())

    abun = rand(Multinomial(individuals, numSpecies))
    trts = TraitCollection2(temp_traits, prec_traits)

    death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
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

    # Load ERA climates
    dir1 = "data/ERA"
    tempax2 = readERA(dir1, "era_int_temp2m_", "t2m", times)
    precax2 = readERA(dir1, "era_int_totalprec_", "tp", times)
    soilwaterax2 = readERA(dir1, "era_int_soilwater1", "swvl1", times)
    solarradax2 = readERA(dir1, "era_int_netsolar", "ssr", times)

    # Concatenate and crop to Africa
    temp = cat(tempax1.array[:, :, 1901year .. (1979year-1month)], tempax2.array, dims = 3)
    prec = cat(precax1.array[:, :, 1901year .. (1979year-1month)], precax2.array, dims = 3)
    prec[prec .< 0m] *= 0
    soilwater = cat(soilwaterax1.array[:, :, 1901year .. (1979year-1month)], soilwaterax2.array, dims = 3)
    soilwater[soilwater .< 0m^3] *= 0
    solarrad = cat(solarradax1.array[:, :, 1901year .. (1979year-1month)], solarradax2.array, dims = 3)
    solarrad[solarrad .< 0J/m^2] *= 0

    africa_temp = ERA(temp[-25°.. 50°, -35°.. 40°, :])
    africa_prec = ERA(prec[-25°.. 50°, -35°.. 40°, :])
    africa_prec.array = AxisArray(uconvert.(mm, africa_prec.array), africa_prec.array.axes)

    # Convert water and solar to the right area/volume for one grid square
    africa_water = uconvert.(m^3, soilwater[-25°.. 50°, -35°.. 40°, :] .* (80km * 80km * 7cm) ./ m^3)
    africa_solar = uconvert.(kJ, solarrad[-25°.. 50°, -35°.. 40°, :] .* (80km * 80km))

    africa_solar = SolarTimeBudget(africa_solar, 1)
    africa_water = EcoSISTEM.VolWaterTimeBudget(africa_water, 1)
    bud = BudgetCollection2(africa_solar, africa_water)
    active =  Array{Bool, 2}(.!isnan.(africa))

    # Put together abiotic environment
    ae1 = eraAE(africa_temp, africa_solar, active)
    ae2 = eraAE(africa_prec, africa_solar, active)
    hab = HabitatCollection2(ae1.habitat, ae2.habitat)
    ae = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, ae1.active, bud, ae1.names)

    rel1 = Gauss{eltype(ae.habitat.h1)}()
    rel2 = Gauss{eltype(ae.habitat.h2)}()
    rel = multiplicativeTR2(rel1, rel2)
    eco = Ecosystem(ecopopulate!, sppl, ae, rel)

    Chapter5.reverseTime!(eco)
    simulate!(eco, 118years - 1month, 1month)
    Chapter5.forwardTime!(eco)
    runburnin!(eco, 10years, 1month)
    return eco
end

output = runSim()
JLD.save("Africa_rand_1.jld", "diver", EcoSISTEM.SavedLandscape(output.abundances))

function buildEco(outputfile::String, repeatYear::Bool, cacheFolder::String, cacheFile::String)
    gbif = JuliaDB.load("data/Full_GBIF_africa")
    traits = JuliaDB.load("data/Africa_traits")
    traits = filter(tr -> (tr.swvl1 > 0.0m^3) & (tr.ssr > 0.0J/m^2) & (tr.prec_sd > 0.0mm) & (tr.tmin_sd > 0.0K), traits)
    gbif = filter(gb -> gb.SppID in select(traits, :SppID), gbif)

    file = "data/Africa.tif"
    africa = readfile(file, -25°, 50°, -35°, 40°)[:, end:-1:1]
    #heatmap(africa)

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
    temp_traits = GaussTrait(tmean, tsd)

    pmean = collect(select(traits, :prec_mean))
    psd = collect(select(traits, :prec_sd))
    prec_traits = GaussTrait(pmean, psd)

    av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
    kernel = GaussianKernel(av_dist, 10e-10)
    movement = BirthOnlyMovement(kernel, NoBoundary())

    abun = rand(Multinomial(individuals, numSpecies))
    trts = TraitCollection2(temp_traits, prec_traits)

    death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
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

    # Load ERA climates
    dir1 = "data/ERA"
    tempax2 = readERA(dir1, "era_int_temp2m_", "t2m", times)
    precax2 = readERA(dir1, "era_int_totalprec_", "tp", times)
    soilwaterax2 = readERA(dir1, "era_int_soilwater1", "swvl1", times)
    solarradax2 = readERA(dir1, "era_int_netsolar", "ssr", times)

    # Concatenate and crop to Africa
    temp = cat(tempax1.array[:, :, 1901year .. (1979year-1month)], tempax2.array, dims = 3)
    prec = cat(precax1.array[:, :, 1901year .. (1979year-1month)], precax2.array, dims = 3)
    prec[prec .< 0m] *= 0
    soilwater = cat(soilwaterax1.array[:, :, 1901year .. (1979year-1month)], soilwaterax2.array, dims = 3)
    soilwater[soilwater .< 0m^3] *= 0
    solarrad = cat(solarradax1.array[:, :, 1901year .. (1979year-1month)], solarradax2.array, dims = 3)
    solarrad[solarrad .< 0J/m^2] *= 0

    africa_temp = ERA(temp[-25°.. 50°, -35°.. 40°, :])
    africa_prec = ERA(prec[-25°.. 50°, -35°.. 40°, :])
    africa_prec.array = AxisArray(uconvert.(mm, africa_prec.array), africa_prec.array.axes)

    # Convert water and solar to the right area/volume for one grid square
    africa_water = uconvert.(m^3, soilwater[-25°.. 50°, -35°.. 40°, :] .* (80km * 80km * 7cm) ./ m^3)
    africa_solar = uconvert.(kJ, solarrad[-25°.. 50°, -35°.. 40°, :] .* (80km * 80km))

    africa_solar = SolarTimeBudget(africa_solar, 1)
    africa_water = EcoSISTEM.VolWaterTimeBudget(africa_water, 1)
    bud = BudgetCollection2(africa_solar, africa_water)
    active =  Array{Bool, 2}(.!isnan.(africa))

    # Put together abiotic environment
    ae1 = eraAE(africa_temp, africa_solar, active)
    ae2 = eraAE(africa_prec, africa_solar, active)
    hab = HabitatCollection2(ae1.habitat, ae2.habitat)
    ae = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, ae1.active, bud, ae1.names)

    rel1 = Gauss{eltype(ae.habitat.h1)}()
    rel2 = Gauss{eltype(ae.habitat.h2)}()
    rel = multiplicativeTR2(rel1, rel2)
    eco = Ecosystem(emptypopulate!, sppl, ae, rel)

    output = JLD.load(outputfile, "diver")
    eco.abundances = GridLandscape(output, (numSpecies, grd[1], grd[2]))
    divfuns = [norm_sub_alpha, raw_sub_alpha, norm_sub_beta, raw_sub_beta, norm_sub_rho, raw_sub_rho, sub_gamma]
    q = [0.0, 1.0, 2.0, Inf]
    simDict = Dict("times" => 118years - 1month, "burnin" => 0year, "interval" => 1year, "timestep" => 1month, "divfuns" => divfuns, "q" => q, "cacheInterval" => 10years, "fileName" => cacheFile)
    lensim = length(0month:simDict["interval"]:simDict["times"])
    diver = zeros(length(divfuns), size(eco.abundances.matrix, 2), lensim, length(q))

    #runburnin!(eco, simDict["burnin"], simDict["timestep"])
    if repeatYear
        Chapter5.keepYear!(eco)
    end
    runsim!(diver, eco, simDict, cacheFolder)
    return diver
end

#output = buildEco("Africa_diversity_1.jld", false, "/media/storage/Chapter5/")
#JLD.save("Africa_full_1.jld", "diver", output)

#output = buildEco("Africa_diversity_1.jld", true, "/media/storage/Chapter5/neutral/")
#JLD.save("Africa_neutral_1.jld", "diver", output)


# ws2
for i in 1:10
    output = buildEco("Africa_rand_1.jld", false, "/media/storage/Chapter5/Africa/rand/run/", "Africa_$i")
    JLD.save("Africa_rand_full_$i.jld", "diver", output)

    output = buildEco("Africa_rand_1.jld", true, "/media/storage/Chapter5/Africa/rand/neutral/", "Africa_$i")
    JLD.save("Africa_rand_neutral_$i.jld", "diver", output)
    print("$i \n")
end
