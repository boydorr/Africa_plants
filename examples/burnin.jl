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
using StatsBase
using EcoSISTEM
using Distributions
using AxisArrays
using Random

function runSim()
    gbif = JuliaDB.load("GBIF_africa_fil")
    gbif = filter(g -> !ismissing(g.date), gbif)
    traits = JuliaDB.load("Africa_traits_fil")
    file = "Africa.tif"
    africa = readfile(file, -25°, 50°, -35°, 40°)[:, end:-1:1]
    #heatmap(africa)

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

    ## Load CERA climates
    dir1 = "CERA"
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
    dir1 = "ERA"
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
    eco = Ecosystem(sppl, ae, rel)

    file = "Africa_effort.tif"
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
JLD.save("1901.jld", "diver", EcoSISTEM.SavedLandscape(output.abundances))
sum(mapslices(sum, output.abundances.matrix, dims =2)[:,1].>0)
