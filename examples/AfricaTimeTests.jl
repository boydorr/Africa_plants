#include("GIT/Chapter5/src/Chapter5.jl")
include("../src/Chapter5.jl")
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
using TimerOutputs
#using Plots
#pyplot()
numSpp =100
#function runSim(numSpp::Int64)
    gbif = JuliaDB.load("data/Full_GBIF_africa")
    traits = JuliaDB.load("data/Full_africa_traits")
    species = collect(JuliaDB.select(traits, :SppID))[1:numSpp]
    traits = filter(tr -> tr.SppID in species, traits)
    gbif = filter(gb -> gb.SppID in species, gbif)
    traits = filter(tr -> !isnan(tr.swvl1) & (tr.swvl1 > 0.0m^3), traits)

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

    trng = collect(select(traits, :trng))
    tmean = collect(select(traits, :tmean))
    temp_traits = GaussTrait(tmean, trng)

    prec = uconvert.(mm, collect(select(traits, :tp)))
    prec_traits = GaussTrait(prec, fill(mean(prec), length(prec)))

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
    soilwater = cat(soilwaterax1.array[:, :, 1901year .. (1979year-1month)], soilwaterax2.array, dims = 3)
    solarrad = cat(solarradax1.array[:, :, 1901year .. (1979year-1month)], solarradax2.array, dims = 3)

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

    start = startingArray(gbif, numSpecies)
    #startabun = reshape(mapslices(sum, start, dims = 1), 100, 100)
    #heatmap(log.(1 .+ startabun))

    adjust1 = [eco.abenv.budget.b1.matrix[x,y,end]/ sum(eco.spplist.requirement.r1.energy .* start[:, convert_coords(x,y, 100)]) for y in 1:100 for x in 1:100]
    adjust1[.!active[1:end]] .= 0

    adjust2 = [eco.abenv.budget.b2.matrix[x,y,end]/sum(eco.spplist.requirement.r2.energy .* start[:, convert_coords(x,y, 100)])  for y in 1:100 for x in 1:100]
    adjust2[.!active[1:end]] .= 0

    adjust = [minimum([adjust1[i], adjust2[i]]) for i in eachindex(adjust1)]
    adjust[isnan.(adjust)] .= 1
    adjust[isinf.(adjust)] .= 1
    #heatmap(reshape(log.(1 .+ adjust), 100, 100))

    eco.abundances.matrix .= round.(Int64, start .* hcat(adjust...))
    #adjustabun = mapslices(sum, eco.abundances.grid, dims = 1)[1, :, :]
    #heatmap(rotl90(log.(1 .+ adjustabun)[:, end:-1:1]))
    eco.abenv.budget.b2.matrix .+= 1e-9m^3
    file = "data/Africa_effort.tif"
    effort = Array(readfile(file, -25°, 50°, -35°, 40°))[:, end:-1:1]
    effort[isnan.(effort)] .= 1

    #fillarray = Array{Int64, 2}(undef, size(eco.abundances.matrix, 1), size(eco.abundances.matrix, 2))
    #heatmap(effort)
    #print("Start sim")
    #simtime = TimerOutput()
    #@timeit simtime "simulate" reverseBurnin!(eco, 1year, 1month, gbif, effort, fillarray)
    #return simtime
#end

output = runSim(100)
#print(Float64(output))


divfuns = [raw_sub_alpha, norm_sub_beta, raw_sub_beta, norm_sub_rho, raw_sub_rho, sub_gamma]
q = 1.0
simDict = Dict("times" => 1years, "burnin" => 0year, "interval" => 1month, "timestep" => 1month, "divfuns" => divfuns, "q" => q, "reps" => 50, "cacheInterval" => 2months, "fileName" => "Africa")
lensim = length(0month:simDict["interval"]:simDict["times"])
diver = zeros(length(divfuns), size(eco.abundances.matrix, 2), lensim)

runsim!(diver, eco, simDict, "data/test")
