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

    # death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
    # birth_rates = death_rates
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
    #effort[effort .> 1e3] .= 1e3

    start = startingArray(gbif, numSpecies, true) .* hcat(effort[1:end]...)

    # start[start .> 0] .= 1

    # adjust1 = [eco.abenv.budget.b1.matrix[x,y,end]/ sum(eco.spplist.requirement.r1.energy .* start[:, convert_coords(x,y, 100)]) for y in 1:100 for x in 1:100]
    # adjust1[.!active[1:end]] .= 0

    # adjust2 = [eco.abenv.budget.b2.matrix[x,y,end]/sum(eco.spplist.requirement.r2.energy .* start[:, convert_coords(x,y, 100)])  for y in 1:100 for x in 1:100]
    # adjust2[.!active[1:end]] .= 0

    # adjust = [minimum([adjust1[i], adjust2[i]]) for i in eachindex(adjust1)]
    # adjust[isnan.(adjust)] .= 1
    # adjust[isinf.(adjust)] .= 1

    # adjust[adjust .< 1] .= 1
    # #adjust[adjust .> 500] .= 500

    # adjusted_start = round.(Int64, start .* hcat(adjust...))
    #adjusted_start[(adjusted_start .< 1e3) .& (adjusted_start .> 0)] .= 1e3

    eco.abundances.matrix .+= round.(Int64, start)
    eco.abenv.budget.b2.matrix .+= 1e-9m^3
    fillarray = Array{Int64, 2}(undef, size(eco.abundances.matrix, 1), size(eco.abundances.matrix, 2))
    reverseBurnin!(eco, 118years - 1month, 1month, 10years, "../../sdb/Chapter5/burnin/", gbif, effort, fillarray)
    runburnin!(eco, 10years, 1month)
    return eco
end

# for i in 1:3
#     output = runSim()
#     JLD.save("Africa_diversity_$i.jld", "diver", EcoSISTEM.SavedLandscape(output.abundances))
# end

output = runSim()
JLD.save("Africa_start2.jld", "diver", EcoSISTEM.SavedLandscape(output.abundances))
sum(mapslices(sum, output.abundances.matrix, dims =2)[:,1].>0)

function buildEco(outputfile::String, repeatYear::Bool, cacheFolder::String, cacheFile::String)
    gbif = JuliaDB.load("GBIF_africa_fil")
    gbif = filter(g -> !ismissing(g.date), gbif)
    traits = JuliaDB.load("Africa_traits_fil")

    file = "Africa.tif"
    africa = readfile(file, -25.0°, 50.0°, -35.0°, 40.0°)[:, end:-1:1]
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
    africa_prec.array .= AxisArray(uconvert.(mm, africa_prec.array), africa_prec.array.axes)

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

#output = buildEco("Africa_diversity_1.jld", false, "/media/storage/Chapter5/")
#JLD.save("Africa_full_1.jld", "diver", output)

#output = buildEco("Africa_diversity_1.jld", true, "/media/storage/Chapter5/neutral/")
#JLD.save("Africa_neutral_1.jld", "diver", output)
output = buildEco("1901.jld", false, "/home/claireh/sdb/Chapter5/run/", "Africa_1")

# ws2
for i in 2:10
    output = buildEco("1901.jld", false, "/home/claireh/sdc/Chapter5/run/", "Africa_$i")
    JLD.save("Africa_new_full_$i.jld", "diver", output)

    output = buildEco("1901.jld", true, "/home/claireh/sdc/Chapter5/neutral/", "Africa_$i")
    JLD.save("Africa_new_neutral_$i.jld", "diver", output)
    print("$i \n")
end

#ws4
for i in 1:10
    output = buildEco("Africa_new_1.jld", false, "Africa/new/run/", "Africa_$i")
    JLD.save("Africa_new_full_$i.jld", "diver", output)

    output = buildEco("Africa_new_1.jld", true, "Africa/new/neutral/", "Africa_$i")
    JLD.save("Africa_new_neutral_$i.jld", "diver", output)
    print("$i \n")
end

#ws3
for i in 1:10
    output = buildEco("Africa_diversity_1.jld", false, "run/", "Africa_$i")
    JLD.save("Africa_full_$i.jld", "diver", output)

    output = buildEco("Africa_diversity_1.jld", true, "neutral/", "Africa_$i")
    JLD.save("Africa_neutral_$i.jld", "diver", output)
    print("$i \n")
end

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.ClimatePref
using JLD
using Plots
using Diversity
using Statistics

function catmean(filename::String, num::Int64)
    diver = JLD.load(filename*"1.jld", "div")
    for i in 1:num
        diveri = JLD.load(filename*"$i.jld", "div")
        diver = cat(diver, diveri, dims = 6)
    end
    return mapslices(mean, diver, dims = 6)[:, :, :, :, :, 1]
end

diver = catmean("data/run/Africa_", 10)
#div = JLD.load("data/Africa_1.jld", "div")
africa = readfile("data/Africa.tif", -25°, 50°, -35°, 40°)[:, end:-1:1]
active =  Array{Bool, 2}(.!isnan.(africa))
function plotProg(div, active, qind, divind; logit=true, kwargs...)
    years = [1, 25, 50, 118]
    h = heatmap(background_color = :lightblue, background_color_outside=:white, grid = false, color = :algae_r, layout = (@layout [a b; c d]), size = (1000, 750), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 5.0*Plots.mm)
    for i in eachindex(years)
        year = div[divind, :, years[i], qind]
        year = reshape(year, size(active))
        year[isnan.(year)] .= 0
        year[.!(active)] .= NaN
        logyear = ifelse(logit, transpose(log.(1 .+ year)), transpose(year))
        h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, logyear, background_color = :lightblue, background_color_outside=:white, grid = false, color = :algae_r, subplot = i; kwargs...)
    end
    display(h)
end
plotProg(diver, active, 1, 1, logit=false)

spr = reshape(diver[1, :, 118, 1] .- diver[1, :, 1, 1], size(active))
spr[isnan.(spr)] .= 0
spr[.!(active)] .= NaN
h = heatmap(-25:0.75:49.25, -35:0.75:39.25, spr', 
background_color = :lightblue, background_color_outside=:white, 
grid = false, color = :algae)

div_neutral = catmean("data/neutral/Africa_", 10)
#div_neutral = JLD.load("data/Africa_1_neutral.jld", "div")
plotProg(div_neutral, active, 1, 1, logit=false)

sr_neutral = reshape(div_neutral[1, :, 118, 1] .- div_neutral[1, :, 1, 1], size(active))
sr_neutral[isnan.(sr_neutral)] .= 0
sr_neutral[.!(active)] .= NaN
h = heatmap(-25:0.75:49.25, -35:0.75:39.25, sr_neutral', 
background_color = :lightblue, background_color_outside=:white, 
grid = false, color = :algae)


h = heatmap(-25:0.75:49.25, -35:0.75:39.25, sr_neutral' .- spr', 
background_color = :lightblue, background_color_outside=:white, 
grid = false, color = cgrad([:darkorange, :white, :darkblue], [0.0, 0.25, 1.0]), clim = (-500, 1500))

function plotDiff(div1, div2, active, qind, divind; kwargs...)
    years = [1, 25, 50, 118]
    div = div1 .- div2
    h = heatmap(background_color = :lightblue, background_color_outside=:white, grid = false, color = :algae_r, layout = (@layout [a b; c d]), size = (1000, 750), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 5.0*Plots.mm)
    for i in eachindex(years)
        year = div[divind, :, years[i], qind]
        year = reshape(year, size(active))
        year[isnan.(year)] .= 0
        year[.!(active)] .= NaN
        h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, year', 
        background_color = :lightblue, background_color_outside=:white, 
        grid = false, color = :algae_r, subplot = i; kwargs...)
    end
    display(h)
end
plotDiff(diver, div_neutral,  active, 1, 1, clim = (-500, 500), color = :balance)


diffs = diver .- div_neutral
diffs = reshape(diffs[1, :, 118], size(active))
diffs[isnan.(diffs)] .= 0
diffs[.!(active)] .= NaN


years = [1, 50, 118]
titles = ["A", "B", "C"]
h = heatmap(background_color = :lightblue, background_color_outside=:white, grid = false, color = :algae_r, layout = (@layout [a b; c d]), size = (1000, 750), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 5.0*Plots.mm)
for i in eachindex(years)
    year = diver[1, :, years[i], 1]
    year = reshape(year, size(active))
    year[isnan.(year)] .= 0
    year[.!(active)] .= NaN
    h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, year', 
    background_color = :lightblue, 
    background_color_outside=:white, 
    title = titles[i], titleloc = :left,
    grid = false, color = :algae_r, subplot = i)
end
display(h)
h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, diffs', 
background_color = :lightblue, background_color_outside=:white, 
grid = false, color = cgrad([:darkblue, :white, :darkorange], [0.0, 0.5, 1.0]),
 clim = (-300, 300), subplot = 4, title = "D", titleloc = :left)
Plots.pdf("plots/AfricaChange.pdf")

years = [1, 118]
titles = ["A", "B"]
h = heatmap(background_color = :lightblue, background_color_outside=:white,
grid = false, color = :algae_r, layout = (@layout [a b]), size = (1000, 600), 
guidefontsize = 16,tickfontsize= 16, titlefontsize=18, 
margin = 5.0*Plots.mm)
for i in eachindex(years)
    year = diver[1, :, years[i], 1]
    year = reshape(year, size(active))
    year[isnan.(year)] .= 0
    year[.!(active)] .= NaN
    h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, year', 
    background_color = :lightblue, 
    background_color_outside=:white, 
    title = titles[i], titleloc = :left,
    grid = false, color = :algae_r, subplot = i, aspect_ratio = 1)
end
display(h)
Plots.pdf("plots/AfricaChangeNew.pdf")


function SpeciesLoss(file1::String)
    abun1 = JLD.load(file1, "div")
    sumabun1 = abun1[1, :, 1, 1]
    sumabun2 = abun1[1, :, 118, 1]
return [sum(sumabun1 .> 0), sum(sumabun2 .> 0)]
end


function catmean(file1)
    spp = zeros(10, 2)
    for j in 1:10
        spp[j, :] = SpeciesLoss(file1*"$j.jld")
        print(j, "\n")
    end
    return spp
end

div_full = catmean("data/run/Africa_")
div_neutral = catmean("data/neutral/Africa_", 10)

diff_tot = div_full[1, :, 118, 1] .- div_full[1, :, 1, 1] 
mean(diff_tot[.!isnan.(diff_tot)])

diff_tot = div_neutral[1, :, 118, 1] .- div_neutral[1, :, 1, 1] 
mean(diff_tot[.!isnan.(diff_tot)])