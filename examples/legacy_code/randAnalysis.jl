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

divfuns = [norm_sub_alpha, raw_sub_alpha, norm_sub_beta, raw_sub_beta, norm_sub_rho, raw_sub_rho, sub_gamma]
q = [0.0, 1.0, 2.0, Inf]
africa = readfile("data/Africa.tif", -25°, 50°, -35°, 40°)[:, end:-1:1]
active =  Array{Bool, 2}(.!isnan.(africa))
function plotProg(div, active, qind, divind; logit=true, kwargs...)
    years = [1, 25, 50, 118]
    div1 = JLD.load("data/Africa_rand_full_1.jld", "diver")
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
        div1 = JLD.load("data/Africa_rand_full_1.jld", "diver")
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
            yr = years[i]
            mar = ifelse(i == 4, 10.0*Plots.mm, 1.0*Plots.mm,)
            title = ifelse(j == 1, "Year $yr", "")
            h = heatmap!(-25:0.75:49.25, -35:0.75:39.25, logyear, background_color = :lightblue, background_color_outside=:white, grid = false, color = :algae_r, subplot = count, clim = (mins, maxs), colorbar = false, right_margin = mar, title = title; kwargs...)
        end
    end
    display(h)
end

diver = JLD.load("data/Africa_new_full_1.jld", "diver")
diver = reshape(diver, 7, 100, 100, 118, 4)

plotProg(diver, active, 2, 1, colorbar_title = "ᾱ", logit = true)
Plots.pdf("plots/rand/RandNormAlphas.pdf")

plotProg(diver, active, 1, 1, colorbar_title = "ᾱ", logit = false, color = cgrad(:algae_r, scale = :exp))
Plots.pdf("plots/rand/RandSpeciesRichness.pdf")

plotProg(diver, active, 2, 3, colorbar_title = "β̄", clim = (0, 10))
Plots.pdf("plots/rand/RandNormBetas.pdf")

plotProg(diver, active, 2, 5, colorbar_title = "ρ̄", logit = false)
Plots.pdf("plots/rand/RandNormRhos.pdf")

plotProg(diver, active, 2, 7, colorbar_title = "γ", logit = false)
Plots.pdf("plots/rand/RandGammas.pdf")
