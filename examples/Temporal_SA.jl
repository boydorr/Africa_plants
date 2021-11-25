using JLD
using Diversity
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using Statistics
using DataFrames
using Random
using Printf

function meanBeta(i::Int64, file1::String, file2::String)
    output1 = JLD.load(file1 * (@sprintf "%02d.jld" i), "abun")
    output2 = JLD.load(file2 * (@sprintf "%02d.jld" i), "abun")

    nonzeros = findall((output1 .> 0) .& (output2 .> 0))
    sites = unique([nonzeros[i][2] for i in 1:length(nonzeros)])
    betas = zeros(length(sites))
    Threads.@threads for j in eachindex(sites)
        betas[j] = norm_meta_beta(Metacommunity([output1[:, sites[j]] output2[:, sites[j]]]), 1.0)[!, :diversity][1]
    end
    return prod(betas)^(1/length(sites))
end
function compbetas()
    betas = zeros(10, 12)
    for j in 1:10
        betas[j, :] = [meanBeta(i, "neutral/cache/SA_$j", "run/cache/SA_$j") for i in 0:11]
        print(j, "\n")
    end
    return betas
end
betas = compbetas()
save("Betas_comp_SA.jld", "betas", betas)

function neutralbetas()
    betas = zeros(10, 12)
    choice = [1:10 (1:10).+1]
    choice[10, 2] = 1
    for j in 1:10
        x = choice[j, 1]; y= choice[j, 2]
        betas[j, :] = [meanBeta(i, "neutral/cache/SA_$x", "neutral/cache/SA_$y") for i in 0:11]
        print(j, "\n")
    end
    return betas
end
betas2 = neutralbetas()
save("Betas_neutral_SA.jld", "betas", betas2)

function fullbetas()
    betas = zeros(10, 12)
    choice = [1:10 (1:10).+1]
    choice[10, 2] = 1
    for j in 1:10
        x = choice[j, 1]; y= choice[j, 2]
        betas[j, :] = [meanBeta(i, "run/cache/SA_$x", "run/cache/SA_$y") for i in 0:11]
        print(j, "\n")
    end
    return betas
end
betas3 = fullbetas()
save("Betas_full_SA.jld", "betas", betas3)

function meanBeta(i::Int64, file1::String, file2::String)
    output1 = JLD.load(file1, "abun")
    output2 = JLD.load(file2 * (@sprintf "%02d.jld" i), "abun")

    nonzeros = findall((output1 .> 0) .& (output2 .> 0))
    sites = unique([nonzeros[i][2] for i in 1:length(nonzeros)])
    betas = zeros(length(sites))
    Threads.@threads for k in eachindex(sites)
        betas[k] = norm_meta_beta(Metacommunity([output1[:, sites[k]] output2[:, sites[k]]]), 1.0)[!, :diversity][1]
    end
    return prod(betas)^(1/length(sites))
end

function neutralbetas()
    betas = zeros(10, 12)
    for j in 1:10
        betas[j, :] = [meanBeta(i, "neutral/cache/SA_$j"*"00.jld", "neutral/cache/SA_$j") for i in 0:11]
        print(j, "\n")
    end
    return betas
end

betas2 = neutralbetas()
save("Betas_neutral_SA_1901.jld", "betas", betas2)

function fullbetas()
    betas = zeros(10, 12)
    for j in 1:10
        betas[j, :] = [meanBeta(i, "run/cache/SA_$j"*"00.jld", "run/cache/SA_$j") for i in 0:11]
        print(j, "\n")
    end
    return betas
end

betas3 = fullbetas()
save("Betas_full_SA_1901.jld", "betas", betas3)


betas = JLD.load("data/Betas_comp_SA.jld", "betas")
betas2 = JLD.load("data/Betas_neutral_SA.jld", "betas")
betas[1] = betas2[1] = 1.001
#betas3 = JLD.load("data/Betas_full_SA.jld", "betas")
plot(0:11, mapslices(mean, betas, dims = 1)[1,:], xlab = "Decade", ylab = "Average number of distinct subcommunities", grid = false, ribbon = 1.96 .* mapslices(std, betas, dims = 1)[1,:]./sqrt(10), fillalpha=.5,guidefontsize =14,tickfontsize=12, legendfontsize=14, label = "Comparison", size = (800, 600))
plot!(0:11, mapslices(mean, betas2[1:9, :], dims = 1)[1,:], ribbon = 1.96 .* mapslices(std, betas2[1:9, :], dims = 1)[1,:]./sqrt(9), fillalpha=.5, guidefontsize =14,tickfontsize=12, legendfontsize=14, margin = 1.0*Plots.mm, size = (800, 600), label = "Neutral model")
#plot!(0:11, mapslices(mean, betas3, dims = 1)[1,:], ribbon = 1.96 .* mapslices(std, betas3, dims = 1)[1,:]./sqrt(15), fillalpha=.5)

Plots.pdf("plots/betaSA.pdf")

betas = JLD.load("data/Betas_neutral_SA_1901.jld", "betas")
betas2 = JLD.load("data/Betas_full_SA_1901.jld", "betas")
plot(1901:10:2018, mapslices(mean, betas, dims = 1)[1,:], xlab = "Decade", ylab = "Average number of distinct subcommunities", grid = false, ribbon = 1.96 .* mapslices(std, betas, dims = 1)[1,:]./sqrt(10), fillalpha=.5, guidefontsize =14,tickfontsize=12, legendfontsize=14, label = "Neutral climate", size = (800, 600), ylim = (1, 1.1))
plot!(1901:10:2018, mapslices(mean, betas2, dims = 1)[1,:], ribbon = 1.96 .* mapslices(std, betas2, dims = 1)[1,:]./sqrt(10), fillalpha=.5, label = "Full climate")
#plot!(0:11, mapslices(mean, betas3, dims = 1)[1,:], ribbon = 1.96 .* mapslices(std, betas3, dims = 1)[1,:]./sqrt(15), fillalpha=.5)

Plots.pdf("plots/betaSA.pdf")


function meanGammas(file::String)
    gammas = zeros(12)
    for i in 0:11
        output = JLD.load(file * (@sprintf "%02d.jld" i), "abun")
        gammas[i+1] = meta_gamma(Metacommunity(output), 1.0)[!, :diversity][1]
        print(i, "\n")
    end
    return gammas
end
gammas1 = [meanGammas("run/cache/SA_$f") for f in 1:10]
gammas2 = [meanGammas("neutral/cache/SA_$f") for f in 1:10]

JLD.save("Gamma_run_SA.jld", "gammas", gammas1)
JLD.save("Gamma_neutral_SA.jld", "gammas", gammas2)

function meanGammas(file::String)
    gammas = zeros(12)
    for i in 0:11
        output = JLD.load(file * (@sprintf "%02d.jld" i), "abun")
        gammas[i+1] = meta_gamma(Metacommunity(output), 0.0)[!, :diversity][1]
        print(i, "\n")
    end
    return gammas
end
gammas1 = [meanGammas("run/cache/SA_$f") for f in 1:10]
gammas2 = [meanGammas("neutral/cache/SA_$f") for f in 1:10]

JLD.save("SR_run_SA.jld", "gammas", gammas1)
JLD.save("SR_neutral_SA.jld", "gammas", gammas2)

gammas = hcat(JLD.load("data/Gamma_neutral_SA.jld", "gammas")...)
gammas2 = hcat(JLD.load("data/Gamma_run_SA.jld", "gammas")...)
plot(1901:10:2018, mapslices(mean, gammas, dims = 2)[:,1], xlab = "Decade", ylab = "Gamma diversity, G", grid = false, ribbon = 1.96 .* mapslices(std, gammas, dims = 2)[:,1]./sqrt(10), fillalpha=.5, guidefontsize =14,tickfontsize=12, legendfontsize=14, label = "Neutral climate", size = (800, 600), ylim = (0, 1400))
plot!(1901:10:2018, mapslices(mean, gammas2, dims = 2)[:,1], ribbon = 1.96 .* mapslices(std, gammas2, dims = 2)[:,1]./sqrt(10), fillalpha=.5, label = "Full climate")
#plot!(0:11, mapslices(mean, betas3, dims = 1)[1,:], ribbon = 1.96 .* mapslices(std, betas3, dims = 1)[1,:]./sqrt(15), fillalpha=.5)

Plots.pdf("plots/GammaSA.pdf")


sr = hcat(JLD.load("data/SR_neutral_SA.jld", "gammas")...)
sr2 = hcat(JLD.load("data/SR_run_SA.jld", "gammas")...)
plot(1901:10:2018, mapslices(mean, sr, dims = 2)[:,1], xlab = "Decade", ylab = "Gamma diversity, G", grid = false, ribbon = 1.96 .* mapslices(std, sr, dims = 2)[:,1]./sqrt(10), fillalpha=.5, guidefontsize =14,tickfontsize=12, legendfontsize=14, label = "Neutral climate", size = (800, 600))
plot!(1901:10:2018, mapslices(mean, sr2, dims = 2)[:,1], ribbon = 1.96 .* mapslices(std, sr2, dims = 2)[:,1]./sqrt(10), fillalpha=.5, label = "Full climate")
#plot!(0:11, mapslices(mean, betas3, dims = 1)[1,:], ribbon = 1.96 .* mapslices(std, betas3, dims = 1)[1,:]./sqrt(15), fillalpha=.5)

Plots.pdf("plots/GammaSR_SA.pdf")
