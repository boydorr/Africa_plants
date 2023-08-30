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

function runSim()
    gbif = JuliaDB.load("data/Full_GBIF_SA")
    traits = JuliaDB.load("data/SA_traits")
    traits = filter(tr -> (tr.swvl1 > 0.0m^3) & (tr.ssr > 0.0J/m^2) & (tr.prec_sd > 0.0mm) & (tr.tmin_sd > 0.0K), traits)
    gbif = filter(gb -> gb.SppID in select(traits, :SppID), gbif)

    file = "data/SouthAmerica.tif"
    sa = readfile(file, -85°, -30.25°, -60.25°, 20°)[:, end:-1:1]
    #heatmap(africa)

    # Set up grid
    numSpecies = length(traits); grd = (107, 73); area = 54e6km^2;
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

    sa_temp = ERA(temp[-85°.. -30.25°, -60.25°.. 20°, :])
    sa_prec = ERA(prec[-85°.. -30.25°, -60.25°.. 20°, :])
    sa_prec.array = AxisArray(uconvert.(mm, sa_prec.array), sa_prec.array.axes)

    # Convert water and solar to the right area/volume for one grid square
    sa_water = uconvert.(m^3, soilwater[-85°.. -30.25°, -60.25°.. 20°, :] .* (80km * 80km * 7cm) ./ m^3)
    sa_solar = uconvert.(kJ, solarrad[-85°.. -30.25°, -60.25°.. 20°, :] .* (80km * 80km))

    sa_solar = SolarTimeBudget(sa_solar, 1)
    sa_water = EcoSISTEM.VolWaterTimeBudget(sa_water, 1)
    bud = BudgetCollection2(sa_solar, sa_water)
    active =  Array{Bool, 2}(.!isnan.(sa))

    # Put together abiotic environment
    ae1 = eraAE(sa_temp, sa_solar, active)
    ae2 = eraAE(sa_prec, sa_solar, active)
    hab = HabitatCollection2(ae1.habitat, ae2.habitat)
    ae = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, ae1.active, bud, ae1.names)

    rel1 = Gauss{eltype(ae.habitat.h1)}()
    rel2 = Gauss{eltype(ae.habitat.h2)}()
    rel = multiplicativeTR2(rel1, rel2)
    eco = Ecosystem(emptypopulate!, sppl, ae, rel)

    start = startingArray(gbif, numSpecies, -85°, -30.25°, -60.25°, 20°)

    adjust1 = [eco.abenv.budget.b1.matrix[x,y,end]/ sum(eco.spplist.requirement.r1.energy .* start[:, convert_coords(x,y, grd[2])]) for y in 1:grd[1] for x in 1:grd[2]]
    adjust1[.!active[1:end]] .= 0

    adjust2 = [eco.abenv.budget.b2.matrix[x,y,end]/sum(eco.spplist.requirement.r2.energy .* start[:, convert_coords(x,y, grd[2])])  for y in 1:grd[1] for x in 1:grd[2]]
    adjust2[.!active[1:end]] .= 0

    adjust = [minimum([adjust1[i], adjust2[i]]) for i in eachindex(adjust1)]
    adjust[isnan.(adjust)] .= 1
    adjust[isinf.(adjust)] .= 1

    eco.abundances.matrix .= round.(Int64, start .* hcat(adjust...))
    eco.abenv.budget.b2.matrix .+= 1e-9m^3
    file = "data/SA_effort.tif"
    effort = Array(readfile(file, -85°, -30.25°, -60.25°, 20°))[:, end:-1:1]
    effort[isnan.(effort)] .= 1
    fillarray = Array{Int64, 2}(undef, size(eco.abundances.matrix, 1), size(eco.abundances.matrix, 2))
    reverseBurnin!(eco, 118years - 1month, 1month, 10years, "SouthAmerica/burnin/", gbif, effort, fillarray, -85°, -30.25°, -60.25°, 20°)
    runburnin!(eco, 10years, 1month)
    return eco
end

output = runSim()

JLD.save("SouthAmerica/SA_diversity_1.jld", "diver", EcoSISTEM.SavedLandscape(output.abundances))

function buildEco(outputfile::String, repeatYear::Bool, cacheFolder::String, cacheFile::String)
    gbif = JuliaDB.load("data/Full_GBIF_SA")
    traits = JuliaDB.load("data/SA_traits")
    traits = filter(tr -> (tr.swvl1 > 0.0m^3) & (tr.ssr > 0.0J/m^2) & (tr.prec_sd > 0.0mm) & (tr.tmin_sd > 0.0K), traits)
    gbif = filter(gb -> gb.SppID in select(traits, :SppID), gbif)

    file = "data/SouthAmerica.tif"
    sa = readfile(file, -85°, -30.25°, -60.25°, 20°)[:, end:-1:1]
    #heatmap(africa)

    # Set up grid
    numSpecies = length(traits); grd = (107, 73); area = 64e6km^2;
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

    sa_temp = ERA(temp[-85°.. -30.25°, -60.25°.. 20°, :])
    sa_prec = ERA(prec[-85°.. -30.25°, -60.25°.. 20°, :])
    sa_prec.array = AxisArray(uconvert.(mm, sa_prec.array), sa_prec.array.axes)

    # Convert water and solar to the right area/volume for one grid square
    sa_water = uconvert.(m^3, soilwater[-85°.. -30.25°, -60.25°.. 20°, :] .* (80km * 80km * 7cm) ./ m^3)
    sa_solar = uconvert.(kJ, solarrad[-85°.. -30.25°, -60.25°.. 20°, :] .* (80km * 80km))

    sa_solar = SolarTimeBudget(sa_solar, 1)
    sa_water = EcoSISTEM.VolWaterTimeBudget(sa_water, 1)
    bud = BudgetCollection2(sa_solar, sa_water)
    active =  Array{Bool, 2}(.!isnan.(sa))

    # Put together abiotic environment
    ae1 = eraAE(sa_temp, sa_solar, active)
    ae2 = eraAE(sa_prec, sa_solar, active)
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
    output = buildEco("SA_diversity_1.jld", false, "/media/storage/Chapter5/SouthAmerica/run/", "SA_$i")
    JLD.save("SA_full_$i.jld", "diver", output)

    output = buildEco("SA_diversity_1.jld", true, "/media/storage/Chapter5/SouthAmerica/neutral/", "SA_$i")
    JLD.save("SA_neutral_$i.jld", "diver", output)
    print("$i \n")
end

#ws4
for i in 1:10
    output = buildEco("SouthAmerica/SA_diversity_1.jld", false, "SouthAmerica/run/", "SA_$i")
    JLD.save("SouthAmerica/SA_full_$i.jld", "diver", output)

    output = buildEco("SouthAmerica/SA_diversity_1.jld", true, "SouthAmerica/neutral/", "SA_$i")
    JLD.save("SouthAmerica/SA_neutral_$i.jld", "diver", output)
    print("$i \n")
end
