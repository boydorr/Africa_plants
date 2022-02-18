using JLD
using Plots
using Statistics
using StatsBase
using Printf
using DataFrames
#plotlyjs()

function catabun(file1::String)
    abun1 = JLD.load(file1 * (@sprintf "%02d.jld" 0), "abun")
    sumabun1 = mapslices(sum, abun1, dims = 2)[:, 1]
    abun2 = JLD.load(file1 * (@sprintf "%02d.jld" 11), "abun")
    sumabun2 = mapslices(sum, abun2, dims = 2)[:, 1]
    return hcat(sumabun1, sumabun2)
end

function catmean(file1)
    abun = catabun(file1*"2")
    for j in 3:10
        abun = abun .+ catabun(file1*"$j")
        print(j, "\n")
    end
    meanabun = abun ./ 10
    return DataFrame(SppID = 1:size(meanabun, 1),  startabun = meanabun[:,1], endabun = meanabun[:,2])
end

abun = catmean("run/cache/Africa_")
neutralabun = catmean("neutral/cache/Africa_")

save("AbunFull.jld", "abun", abun)
save("AbunNeutral.jld", "abun", neutralabun)


abun = catmean("run/cache/SA_")
save("/home/claireh/Documents/Chapter5/AbunFullSA.jld", "abun", abun)

neutralabun = catmean("neutral/cache/SA_")
save("/home/claireh/Documents/Chapter5/AbunNeutralSA.jld", "abun", neutralabun)


abun = load("data/AbunFull.jld2", "abun")
neutralabun = load("data/AbunNeutral.jld2", "neutralabun")
abun = abun[abun[:startabun] .!= 0, :]
neutralabun = neutralabun[neutralabun[:startabun] .!= 0, :]

plotlyjs()
abunorder = sortperm(abun, :startabun, rev = true)
bar(abun[abunorder, :startabun], grid = false, xlab = "Ranked species", 
ylab = "Abundance 1901", size = (1200, 800), 
guidefontsize = 12,tickfontsize= 12, titlefontsize=18, 
margin = 2.0*Plots.mm, legendfontsize = 12, subplot = 1, 
left_margin = 5.0*Plots.mm, color = 2, 
layout = (@layout [a; b]), label = "",ylim = (1, 1e12), 
linecolor = :match, guide_position = :top, 
bottom_margin = 0.0*Plots.mm, yscale = :log10)
bar!(abun[abunorder, :endabun], xlab = "Ranked species", 
ylab = "Abundance 2018", ylim = (1, 1e12), subplot = 2, 
grid = false, guidefontsize = 12, tickfontsize= 12, 
titlefontsize=18, margin = 2.0*Plots.mm, 
legendfontsize = 12, color =2, linecolor = :match,
 yflip = true, guide_position = :top, link = :x, 
 top_margin = 0.0*Plots.mm, label = "", yscale = :log10)
plot!(abun[abunorder, :startabun], subplot = 2, color = :black, label = "", linewidth = 0.5, yscale = :log10)
Plots.pdf("plots/ClimateAbun.pdf")

abunorder = sortperm(neutralabun, :startabun, rev = true)
bar(neutralabun[abunorder, :startabun], grid = false, xlab = "Ranked species", ylab = "Abundance 1901", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, subplot = 1, left_margin = 5.0*Plots.mm, color = 1, layout = (@layout [a; b]), label = "",ylim = (1, 1e12), linecolor = :match, guide_position = :top, bottom_margin = 0.0*Plots.mm, yscale = :log10)
bar!(neutralabun[abunorder, :endabun], xlab = "Ranked species", ylab = "Abundance 2018", ylim = (1, 1e12), subplot = 2, grid = false, guidefontsize = 12, tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, color =1, linecolor = :match, yflip = true, guide_position = :top, link = :x, top_margin = 0.0*Plots.mm, label = "", yscale = :log10)
plot!(neutralabun[abunorder, :startabun], subplot = 2, color = :black, label = "", linewidth = 0.5, yscale = :log10)
Plots.pdf("plots/NeutralAbun.pdf")

using StatsPlots
abun = JLD.load("data/AbunFull.jld", "abun")
neutralabun = load("data/AbunNeutral.jld", "abun")
traits = JuliaDB.load("data/Africa_traits")
traits = filter(tr -> (tr.swvl1 > 0.0m^3) & (tr.ssr > 0.0J/m^2) & (tr.prec_sd > 0.0mm) & (tr.tmin_sd > 0.0K), traits)
abun[:temp] = JuliaDB.select(traits, :tmin_mean)
neutralabun[:temp] = JuliaDB.select(traits, :tmin_mean)
abun = abun[abun[:startabun] .!= 0, :]
neutralabun = neutralabun[neutralabun[:startabun] .!= 0, :]
@df neutralabun histogram(ustrip.(uconvert.(°C, :temp)), weights = :endabun, bins = 20, xlab = "Temperature preference (°C)", ylab = "Abundance 2018", grid = false, guidefontsize = 14, tickfontsize= 14, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 14, color =1, top_margin = 0.0*Plots.mm, label = "", ylim = (0, 6e12), size = (800, 1000), layout = (@layout [a; b]), link = :x)
@df abun histogram!(ustrip.(uconvert.(°C, :temp)), weights = :endabun, bins = 20, xlab = "Temperature preference (°C)", ylab = "Abundance 2018", grid = false, guidefontsize = 14, tickfontsize= 14, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 14, color =2, subplot = 2, top_margin = 0.0*Plots.mm, label = "", ylim = (0, 6e12))
Plots.pdf("plots/TempAbun.pdf")

abun = JLD.load("data/AbunFullSA.jld", "abun")
neutralabun = load("data/AbunNeutralSA.jld", "abun")
abun = abun[abun[:startabun] .!= 0, :]
neutralabun = neutralabun[neutralabun[:startabun] .!= 0, :]

abunorder = sortperm(abun, :startabun, rev = true)
samp = abunorder[1:5:end]
bar(abun[samp, :startabun], grid = false, xlab = "Ranked species", ylab = "Abundance 1901", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, subplot = 1, left_margin = 5.0*Plots.mm, color = 2, layout = (@layout [a; b]), label = "", linecolor = :match, guide_position = :top, bottom_margin = 0.0*Plots.mm, yscale = :log10, ylim = (1, 1e12))
bar!(abun[samp, :endabun], xlab = "Ranked species", ylab = "Abundance 2018", subplot = 2, grid = false, guidefontsize = 12, tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, color =2, linecolor = :match, yflip = true, guide_position = :top, link = :x, top_margin = 0.0*Plots.mm, label = "", yscale = :log10, ylim = (1, 1e12))
plot!(abun[samp, :startabun], subplot = 2, color = :black, label = "", linewidth = 0.5, yscale = :log10)
Plots.pdf("plots/ClimateAbunSA.pdf")

abunorder = sortperm(neutralabun, :startabun, rev = true)
samp = abunorder[1:5:end]
bar(neutralabun[samp, :startabun], grid = false, xlab = "Ranked species", ylab = "Abundance 1901", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, subplot = 1, left_margin = 5.0*Plots.mm, color = 1, layout = (@layout [a; b]), label = "",ylim = (1, 1e12), linecolor = :match, guide_position = :top, bottom_margin = 0.0*Plots.mm, yscale = :log10)
bar!(neutralabun[samp, :endabun], xlab = "Ranked species", ylab = "Abundance 2018", ylim = (1, 1e12), subplot = 2, grid = false, guidefontsize = 12, tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, color =1, linecolor = :match, yflip = true, guide_position = :top, link = :x, top_margin = 0.0*Plots.mm, label = "", yscale = :log10)
plot!(neutralabun[samp, :startabun], subplot = 2, color = :black, label = "", linewidth = 0.5, yscale = :log10)
Plots.pdf("plots/NeutralAbunSA.pdf")

abun = JLD.load("data/AbunFullSA.jld", "abun")
neutralabun = load("data/AbunNeutralSA.jld", "abun")
traits = JuliaDB.load("data/SA_traits")
traits = filter(tr -> (tr.swvl1 > 0.0m^3) & (tr.ssr > 0.0J/m^2) & (tr.prec_sd > 0.0mm) & (tr.tmin_sd > 0.0K), traits)
abun[:temp] = JuliaDB.select(traits, :tmin_mean)
neutralabun[:temp] = JuliaDB.select(traits, :tmin_mean)
abun = abun[abun[:startabun] .!= 0, :]
neutralabun = neutralabun[neutralabun[:startabun] .!= 0, :]
@df neutralabun histogram(ustrip.(uconvert.(°C, :temp)), weights = :endabun, bins = 20, xlab = "Temperature preference (°C)", ylab = "Abundance 2018", size = (800, 1000), grid = false, guidefontsize = 14, tickfontsize= 14, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 14, color =1, top_margin = 0.0*Plots.mm, label = "", ylim = (0, 2e12), layout = (@layout [a; b]), link = :x)
@df abun histogram!(ustrip.(uconvert.(°C, :temp)), weights = :endabun, bins = 20, xlab = "Temperature preference (°C)", ylab = "Abundance 2018", grid = false, guidefontsize = 14, tickfontsize= 14, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 14, color =2, top_margin = 0.0*Plots.mm, label = "", ylim = (0, 2e12), subplot = 2)
Plots.pdf("plots/TempAbunSA.pdf")



function catsum(file1::String)
    abun = zeros(12)
    for i in 0:11
        abun[i+1] = sum(JLD.load(file1 * (@sprintf "%02d.jld" i), "abun"))
    end
    return abun
end

function catmean(file1)
    abun = catsum(file1*"1")
    for j in 2:10
        abun = abun .+ catsum(file1*"$j")
        print(j, "\n")
    end
    return abun ./ 10
end

abun = catmean("cache/Africa_")
save("TotalAbunFull.jld", "abun", abun)
neutralabun = catmean("neutral/cache/Africa_")
save("AbunNeutral.jld", "abun", neutralabun)

abun = catmean("run/cache/SA_")
save("/home/claireh/Documents/Chapter5/TotalAbunFullSA.jld", "abun", abun)
neutralabun = catmean("neutral/cache/SA_")
save("/home/claireh/Documents/Chapter5/TotalAbunNeutralSA.jld", "abun", abun)

function SpeciesLoss(file1::String)
    abun1 = JLD.load(file1 * (@sprintf "%02d.jld" 0), "abun")
    sumabun1 = mapslices(sum, abun1, dims = 2)[:, 1]
    abun2 = JLD.load(file1 * (@sprintf "%02d.jld" 11), "abun")
    sumabun2 = mapslices(sum, abun2, dims = 2)[:, 1]
return [sum(sumabun1 .> 0), sum(sumabun2 .> 0)]
end


function catmean(file1)
    spp = zeros(10, 2)
    for j in 1:10
        spp[j, :] = SpeciesLoss(file1*"$j")
        print(j, "\n")
    end
    return spp
end

abun = catmean("runcache/Africa_")
save("SRChangeFull.jld", "abun", abun)

neutralabun = catmean("neutral/cache/Africa_")
save("SRChangeNeutral.jld", "abun", neutralabun)


abun = catmean("run/cache/SA_")
save("/home/claireh/Documents/Chapter5/SRChangeFullSA.jld", "abun", abun)

neutralabun = catmean("neutral/cache/SA_")
save("/home/claireh/Documents/Chapter5/SRChangeNeutralSA.jld", "abun", neutralabun)

using Plots
plotlyjs()

gbif = JuliaDB.load("data/Full_GBIF_africa")
traits = JuliaDB.load("data/Africa_traits")
traits = filter(tr -> (tr.swvl1 > 0.0m^3) & (tr.ssr > 0.0J/m^2) & (tr.prec_sd > 0.0mm) & (tr.tmin_sd > 0.0K), traits)
gbif = filter(gb -> gb.SppID in JuliaDB.select(traits, :SppID), gbif)
start = startingArray(gbif, length(traits))
startabun = mapslices(sum, start, dims = 2)[:,1]
bar(sort(startabun, rev = true)[1:5:end], grid = false, xlab = "Ranked species", ylab = "GBIF Abundance", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, left_margin = 5.0*Plots.mm, color = 1, label = "", linecolor = :match, guide_position = :top, bottom_margin = 0.0*Plots.mm, yscale = :log10)
Plots.pdf("plots/AfricaAbun.pdf")

gbif = JuliaDB.load("data/Full_GBIF_SA")
traits = JuliaDB.load("data/SA_traits")
traits = filter(tr -> (tr.swvl1 > 0.0m^3) & (tr.ssr > 0.0J/m^2) & (tr.prec_sd > 0.0mm) & (tr.tmin_sd > 0.0K), traits)
gbif = filter(gb -> gb.SppID in select(traits, :SppID), gbif)
start = startingArray(gbif, length(traits), -85°, -30.25°, -60.25°, 20°)
startabun = mapslices(sum, start, dims = 2)[:,1]
bar(sort(startabun, rev = true)[1:5:end], grid = false, xlab = "Ranked species", ylab = "GBIF Abundance", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, left_margin = 5.0*Plots.mm, color = 1, label = "", linecolor = :match, guide_position = :top, bottom_margin = 0.0*Plots.mm, yscale = :log10)
Plots.pdf("plots/SAmericaAbun.pdf")
