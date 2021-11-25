using JLD
using Diversity
using Plots
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using EcoSISTEM.ClimatePref
using Statistics
using EcoSISTEM
using AxisArrays
using Random
plotlyjs()

output = JLD.load("data/1901.jld", "1901")
africa_abun = norm_sub_alpha(Metacommunity(output), 0.0)[:diversity]
africa_abun = transpose(reshape(africa_abun, 100, 100))
africa = readfile("data/Africa.tif", -25°, 50°, -35°, 40°)[:, end:-1:1]
active =  Array{Bool, 2}(.!isnan.(africa))
africa_abun[isnan.(africa_abun)] .= 0
africa_abun[.!(transpose(active))] .= NaN

heatmap(africa_abun, background_color = :lightblue, background_color_outside=:white, grid = false, color = cgrad(:algae, scale = :exp), aspect_ratio = 1)
Plots.pdf("plots/1901.pdf")

output = EcoSISTEM.GridLandscape(JLD.load("Africa_new_1.jld", "diver"), (41844, 100, 100))
africa_abun = norm_sub_alpha(Metacommunity(output.matrix), 0.0)[:diversity]
africa_abun = transpose(reshape(africa_abun, 100, 100))
africa = readfile("data/Africa.tif", -25°, 50°, -35°, 40°)[:, end:-1:1]
active =  Array{Bool, 2}(.!isnan.(africa))
africa_abun[isnan.(africa_abun)] .= 0
africa_abun[.!(transpose(active))] .= NaN

heatmap(africa_abun, background_color = :lightblue, background_color_outside=:white, grid = false, color = cgrad(:algae, scale = :exp), aspect_ratio = 1)
Plots.pdf("Rand1901.pdf")

##
divfuns = [norm_sub_alpha, raw_sub_alpha, norm_sub_beta, raw_sub_beta, norm_sub_rho, raw_sub_rho, sub_gamma]
q = [0.0, 1.0, 2.0, Inf]
africa = readfile("data/Africa.tif", -25°, 50°, -35°, 40°)[:, end:-1:1]
active =  Array{Bool, 2}(.!isnan.(africa))
function plotProg(div, active, qind, divind; logit=true, kwargs...)
    years = [1, 25, 50, 118]
    div1 = JLD.load("data/Africa_full_1.jld", "diver")
    div1 = reshape(div1, 7, 100, 100, 118, 4)
    inds = div1[divind, :, :, years, qind]
    inds[isnan.(inds)] .= 0
    for i in 1:size(inds, 3)
        ind = @view inds[:,:, i]
        ind[.!(active)] .= NaN
    end
    maxs = ifelse(logit, maximum(log.(1 .+ inds[.!isnan.(inds)])), maximum(inds[.!isnan.(inds)]))
    mins = ifelse(logit, minimum(log.(1 .+ inds[.!isnan.(inds)])), minimum(inds[.!isnan.(inds)]))
    h = heatmap(background_color = :lightblue, background_color_outside=:white, grid = false, color = :algae_r, layout = (@layout [a b; c d]), size = (1000, 750), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 5.0*Plots.mm, clim = (mins, maxs))
    for i in eachindex(years)
        year = div[divind, :, :, years[i], qind]
        year[isnan.(year)] .= 0
        year[.!(active)] .= NaN
        logyear = ifelse(logit, transpose(log.(1 .+ year)), transpose(year))
        h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, logyear, background_color = :lightblue, background_color_outside=:white, grid = false, color = :algae_r, subplot = i; kwargs...)
    end
    display(h)
end


function plotProg(div, active, qind;logit=true, kwargs...)
    l = (@layout [a b c d; e f g h; i j k l; m n o p])
    years = [1, 25, 50, 118]
    h = heatmap(background_color = :lightblue, background_color_outside=:white, grid = false, color = :algae_r, layout = l, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm)
    divs = [1, 3, 5, 7]
    count =0
    for j in eachindex(divs)
        div1 = JLD.load("data/Africa_full_1.jld", "diver")
        div1 = reshape(div1, 7, 100, 100, 118, 4)
        inds = div1[divs[j], :, :, years, qind]
        inds[isnan.(inds)] .= 0
        for k in 1:size(inds, 3)
            ind = @view inds[:,:, k]
            ind[.!(active)] .= NaN
        end
        maxs = ifelse(logit, maximum(log.(1 .+ inds[.!isnan.(inds)])), maximum(inds[.!isnan.(inds)]))
        mins = ifelse(logit, minimum(log.(1 .+ inds[.!isnan.(inds)])), minimum(inds[.!isnan.(inds)]))
        for i in eachindex(years)
            count += 1
            year = div[divs[j], :, :, years[i], qind]
            year[isnan.(year)] .= 0
            year[.!(active)] .= NaN
            logyear = ifelse(logit, transpose(log.(1 .+ year)), transpose(year))
            yr = years[i] + 1900
            mar = ifelse(i == 4, 10.0*Plots.mm, 1.0*Plots.mm,)
            title = ifelse(j == 1, "$yr", "")
            h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, logyear, background_color = :lightblue, background_color_outside=:white, grid = false, color = :algae_r, subplot = count, clim = (mins, maxs), colorbar = false, right_margin = mar, title = title; kwargs...)
        end
    end
    display(h)
end
plotlyjs()

function catmean(filename::String, num::Int64)
    diver = JLD.load(filename*"1.jld", "diver")
    diver = reshape(diver, 7, 100, 100, 118, 4)
    for i in 2:num
        diveri = JLD.load(filename*"$i.jld", "diver")
        diveri = reshape(diveri, 7, 100, 100, 118, 4)
        diver = cat(diver, diveri, dims = 6)
    end
    return mapslices(mean, diver, dims = 6)[:, :, :, :, :, 1]
end
diver = catmean("data/Africa_full_", 10)

plotProg(diver, active, 2, 1, colorbar_title = "ᾱ")
Plots.pdf("plots/NormAlphas.pdf")

plotProg(diver, active, 1, 1, colorbar_title = "ᾱ", logit = false)
Plots.pdf("plots/SpeciesRichness.pdf")

plotProg(diver, active, 2, 3, colorbar_title = "β̄")
Plots.pdf("plots/NormBetas.pdf")

plotProg(diver, active, 2, 5, colorbar_title = "ρ̄", logit = false)
Plots.pdf("plots/NormRhos.pdf")

plotProg(diver, active, 2, 7, colorbar_title = "γ")
Plots.pdf("plots/Gammas.pdf")

pyplot()
plotProg(diver, active, 2, 1, colorbar_title = "ᾱ")
Plots.pdf("plots/CB_alphas.pdf")

plotProg(diver, active, 2, 3, colorbar_title = "β̄")
Plots.pdf("plots/CB_betas.pdf")

plotProg(diver, active, 2, 5, colorbar_title = "ρ̄", logit = false)
Plots.pdf("plots/CB_rhos.pdf")

plotProg(diver, active, 2, 7, colorbar_title = "γ")
Plots.pdf("plots/CB_gammas.pdf")

plotlyjs()
plotProg(diver, active, 2)
Plots.pdf("plots/Total.pdf")

plotlyjs()
diver = catmean("data/Africa_neutral_", 10)

plotProg(diver, active, 2, 1, colorbar_title = "ᾱ")
Plots.pdf("plots/NormAlphas_neutral.pdf")

plotProg(diver, active, 2, 3, colorbar_title = "β̄")
Plots.pdf("plots/NormBetas_neutral.pdf")

plotProg(diver, active, 2, 5, colorbar_title = "ρ̄", logit = false)
Plots.pdf("plots/NormRhos_neutral.pdf")

plotProg(diver, active, 2, 7, colorbar_title = "γ")
Plots.pdf("plots/Gammas_neutral.pdf")

plotProg(diver, active, 2)
Plots.pdf("plots/Total_neutral.pdf")


function catstd(filename::String, num::Int64)
    diver = JLD.load(filename*"1.jld", "diver")
    diver = reshape(diver, 7, 100, 100, 118, 4)
    for i in 2:num
        diveri = JLD.load(filename*"$i.jld", "diver")
        diveri = reshape(diveri, 7, 100, 100, 118, 4)
        diver = cat(diver, diveri, dims = 6)
    end
    return mapslices(std, diver, dims = 6)[:, :, :, :, :, 1]
end
diver1 = catstd("data/Africa_full_", 10)
diver2 = catstd("data/Africa_neutral_", 10)
alphas = diver1[1,:, :,:,2]
histogram(alphas[.!isnan.(alphas)], bins = range(0, stop = 12, length = 50), layout = (@layout [a; b]), ylim = (0, 1.5e5))
alphas = diver2[1,:, :,:,2]
histogram!(alphas[.!isnan.(alphas)], bins = range(0, stop = 12, length = 50), subplot = 2, ylim = (0, 1.5e5), color = :green)

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

pyplot()
tempchange = mapslices(mean, africa_temp.array[:, :, 709:end]./K, dims =3)[:, :, 1] .- mapslices(mean, africa_temp.array[:, :, 1:708]./K, dims =3)[:, :, 1]
tempchange[.!(active)] .= NaN
heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(tempchange), layout = (@layout [a b;c d]), subplot = 1, background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1000, 750), title = "Temperature (K)", clim = (-3, 3))
precchange = mapslices(mean, africa_prec.array[:, :, 709:end]./mm, dims =3)[:, :, 1] .- mapslices(mean, africa_prec.array[:, :, 1:708]./mm, dims =3)[:, :, 1]
precchange[.!(active)] .= NaN
heatmap!(-25:0.75:49.25, -35:0.75:39.25, transpose(precchange), layout = (@layout [a b;c d]), subplot = 2, background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1000, 750), title = "Precipitation (mm)", clim = (-4, 4))
waterchange = mapslices(mean, africa_water[:, :, 709:end]./m^3, dims =3)[:, :, 1] .- mapslices(mean, africa_water[:, :, 1:708]./m^3, dims =3)[:, :, 1]
waterchange[.!(active)] .= NaN
heatmap!(-25:0.75:49.25, -35:0.75:39.25, transpose(waterchange), layout = (@layout [a b;c d]), subplot = 3, background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1000, 750), title = "Soil water volume (m³)", clim = (-7e7, 7e7))
solarchange = mapslices(mean, africa_solar[:, :, 709:end]./kJ, dims =3)[:, :, 1] .- mapslices(mean, africa_solar[:, :, 1:708]./kJ, dims =3)[:, :, 1]
solarchange[.!(active)] .= NaN
heatmap!(-25:0.75:49.25, -35:0.75:39.25, transpose(solarchange), layout = (@layout [a b;c d]), subplot = 4, background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1000, 750), title = "Solar radiation (kJ)", clim = (-2e13, 2e13))
Plots.pdf("plots/ClimateChanges50.pdf")

pyplot()
tempchange = mapslices(mean, africa_temp.array[:, :, 1296:1415]./K, dims =3)[:, :, 1] .- mapslices(mean, africa_temp.array[:, :, 1176:1295]./K, dims =3)[:, :, 1]
tempchange[.!(active)] .= NaN
heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(tempchange), layout = (@layout [a b;c d]), subplot = 1, background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1000, 750), title = "Temperature (K)", clim = (-1, 1))
precchange = mapslices(mean, africa_prec.array[:, :, 1296:1415]./mm, dims =3)[:, :, 1] .- mapslices(mean, africa_prec.array[:, :, 1176:1295]./mm, dims =3)[:, :, 1]
precchange[.!(active)] .= NaN
heatmap!(-25:0.75:49.25, -35:0.75:39.25, transpose(precchange), layout = (@layout [a b;c d]), subplot = 2, background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1000, 750), title = "Precipitation (mm)", clim = (-3, 3))
waterchange = mapslices(mean, africa_water[:, :, 1296:1415]./m^3, dims =3)[:, :, 1] .- mapslices(mean, africa_water[:, :, 1176:1295]./m^3, dims =3)[:, :, 1]
waterchange[.!(active)] .= NaN
heatmap!(-25:0.75:49.25, -35:0.75:39.25, transpose(waterchange), layout = (@layout [a b;c d]), subplot = 3, background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1000, 750), title = "Soil water volume (m³)", clim = (-3e7, 3e7))
solarchange = mapslices(mean, africa_solar[:, :, 1296:1415]./kJ, dims =3)[:, :, 1] .- mapslices(mean, africa_solar[:, :, 1176:1295]./kJ, dims =3)[:, :, 1]
solarchange[.!(active)] .= NaN
heatmap!(-25:0.75:49.25, -35:0.75:39.25, transpose(solarchange), layout = (@layout [a b;c d]), subplot = 4, background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1000, 750), title = "Solar radiation (kJ)", clim = (-8e12, 8e12))
Plots.pdf("plots/ClimateChanges10.pdf")

using Random
using EcoSISTEM
using JLD
using EcoSISTEM.ClimatePref
using Unitful
using Unitful.DefaultSymbols
using Plots
plotlyjs()

for i in 1:3
    output = EcoSISTEM.GridLandscape(JLD.load("Africa_diversity_$i.jld", "diver"), (41844, 100, 100))
    africa = readfile("data/Africa.tif", -25°, 50°, -35°, 40°)[:, end:-1:1]
    active =  Array{Bool, 2}(.!isnan.(africa))
    abun = mapslices(sum, output.grid, dims = 1)[1, :, :] .* 1.0
    abun[isnan.(abun)] .= 0.0
    abun[.!(active)] .= NaN
    if i ==1
        h = heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(abun), background_color = :lightblue, background_color_outside=:white, grid = false, color = :algae_r, layout = (@layout [a b; c d]), size = (1000, 750), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 5.0*Plots.mm, subplot = 1, clim = (0, 3e10))
    else
        h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, transpose(abun), subplot = i, grid = false, color = :algae_r, guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 5.0*Plots.mm, clim = (0, 3e10))
    end
end
Plots.pdf("Burnins.pdf")

function catreps(filename::String, num::Int64)
    diver = JLD.load(filename*"1.jld", "diver")
    diver = reshape(diver, 7, 100, 100, 118, 4)
    for i in 2:num
        diveri = JLD.load(filename*"$i.jld", "diver")
        diveri = reshape(diveri, 7, 100, 100, 118, 4)
        diver = cat(diver, diveri, dims = 6)
    end
    return diver
end
diver = catreps("data/Africa_full_", 10)
diverN = catreps("data/Africa_neutral_", 10)
#slopemat = slopediv(diver[:, :, :, :, 2, :]) .* 10
slopemat = slopediv(diver) .* 10
slopematN = slopediv(diverN) .* 10
meanslope = mapslices(mean, slopemat, dims = 5)[:, :, :, :, 1]
meanslopeN = mapslices(mean, slopematN, dims = 5)[:, :, :, :, 1]
meandiff = meanslope .- meanslopeN
#meanslope[:, 20, 84] .= 0
#meanslope[:, 98, 29] .= 0
pyplot()
for i in 1:4
    normalpha = meandiff[1, :, :, i]
    normalpha[isnan.(normalpha)] .= 0
    normalpha[.!(active)] .= NaN
    qs = [0,1,2,Inf]
    maxs = maximum(normalpha[.!isnan.(normalpha)])
    mins = minimum(normalpha[.!isnan.(normalpha)])
    lims = ifelse(abs(maxs) > abs(mins), (-maxs, maxs), (mins, -mins))
    if i ==1
        q = Int64(qs[i])
        h = heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(normalpha),background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm,clim = lims, layout = (@layout [a b;c d]), title = "q = $q")
    else
        if i < 4
            q = Int64(qs[i])
        else
            q = Inf
        end
        title = ifelse(i==4, "q = ∞", "q = $q")
        h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, transpose(normalpha),background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm, clim = lims, subplot = i, title = title)
    end
    display(h)
end
Plots.png("plots/ChangeNormAlphas.png")
normalpha = meanslope[1, :, :]
normalpha[isnan.(normalpha)] .= 0
normalpha[.!(active)] .= NaN
h = heatmap(transpose(normalpha),background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm, clim = (-4, 4))

normbeta = meanslope[3, :, :, 2]
normbeta[isnan.(normbeta)] .= 0
normbeta[.!(active)] .= NaN
h = heatmap(transpose(normbeta),background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm, clim = (-7e4, 7e4))

beta = meanslope[4, :, :]
beta[isnan.(beta)] .= 0
beta[.!(active)] .= NaN
h = heatmap(transpose(beta),background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm, clim = (-0.05, 0.05))

normrho = meanslope[5, :, :, 2]
normrho[isnan.(normrho)] .= 0
normrho[.!(active)] .= NaN
h = heatmap(transpose(normrho),background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm, clim = (-0.005, 0.005))

for i in 1:4
    normrho = meandiff[5, :, :, i]
    normrho[isnan.(normrho)] .= 0
    normrho[.!(active)] .= NaN
    qs = [0,1,2,Inf]
    maxs = maximum(normrho[.!isnan.(normrho)])
    mins = minimum(normrho[.!isnan.(normrho)])
    lims = ifelse(abs(maxs) > abs(mins), (-maxs, maxs), (mins, -mins))
    if i ==1
        q = Int64(qs[i])
        h = heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(normrho),background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm, clim = lims, layout = (@layout [a b;c d]), title = "q = $q",minorgrid=false)
    else
        if i < 4
            q = Int64(qs[i])
        else
            q = Inf
        end
        title = ifelse(i==4, "q = ∞", "q = $q")
        h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, transpose(normrho),background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm, clim = lims, subplot = i, title = title)
    end
    display(h)
end
Plots.png("plots/ChangeNormRhos.png")

rho = meanslope[6, :, :]
rho[isnan.(rho)] .= 0
rho[.!(active)] .= NaN
h = heatmap(transpose(rho),background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm, clim = (-5e2, 5e2))

for i in 1:100
    for j in 1:100
        gammacat[i, j] = floor()
    end
end
gamma = meanslope[7, :, :]
gamma[isnan.(gamma)] .= 0
gamma[.!(active)] .= NaN
h = heatmap(transpose(gamma),background_color = :lightblue, background_color_outside=:white, grid = false, color = :RdBu, size = (1500, 1125), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 1.0*Plots.mm, clim = (-1e3, 1e3))

meanslope = mapslices(mean, slopemat, dims = 4)[:, :, :, 1]
