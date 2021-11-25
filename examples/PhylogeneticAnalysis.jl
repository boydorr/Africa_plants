using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using JuliaDB
using Diversity
using PhyloNetworks
using JLD
using DataFrames

traits = JuliaDB.load("data/Africa_traits")
traits = filter(tr -> (tr.swvl1 > 0.0m^3) & (tr.ssr > 0.0J/m^2) & (tr.prec_sd > 0.0mm) & (tr.tmin_sd > 0.0K), traits)
tree = readTopology("../Chapter4/Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")

SppDict = JLD.load("data/SppDict.jld", "SppDict")
SppID = Dict(value => key for (key, value) in SppDict)
spp_names = [SppDict[x] for x in JuliaDB.select(traits, :SppID)]

cross_species = spp_names ∩ tip_names
phylo_ids = [SppID[x] for x in cross_species]
trait_ids = JuliaDB.select(traits, :SppID)
new_phylo_ids = indexin(phylo_ids, trait_ids)
cross_species = cross_species[.!isnothing.(new_phylo_ids)]
new_phylo_ids = new_phylo_ids[.!isnothing.(new_phylo_ids)]

# Load tree data and select names of species
t = open(io -> parsenewick(io, NamedPolytomousTree),"/home/claireh/Documents/Chapter4/Qian2016.tree")
keeptips!(t, join.(split.(cross_species, " "), "_"))
ty = PhyloTypes(t)

function meanBeta(i::Int64, file1::String, file2::String)
    output1 = JLD.load(file1 * (@sprintf "%02d.jld" i), "abun")
    output2 = JLD.load(file2 * (@sprintf "%02d.jld" i), "abun")

    nonzeros = findall((output1 .> 0) .& (output2 .> 0))
    sites = unique([nonzeros[i][2] for i in 1:length(nonzeros)])
    betas = 1
    for j in eachindex(sites)[1:10]
        betas *= norm_meta_beta(Metacommunity([output1[new_phylo_ids, sites[j]] output2[new_phylo_ids, sites[j]]], t), 1.0)[!, :diversity][1]
    end
    return betas^(1/length(sites))
end

betas = zeros(6, 12)
for j in 1:6
    betas[j, :] = [meanBeta(i, "neutral/cache/Africa_$j", "cache/Africa_$j") for i in 0:11]
end

file1 =  "/media/storage/Chapter5/Africa/neutral/cache/Africa_1"
file2 =  "/media/storage/Chapter5/Africa/cache/Africa_1"
output1 = JLD.load(file1 * (@sprintf "%02d.jld" 1), "abun")
output2 = JLD.load(file2 * (@sprintf "%02d.jld" 1), "abun")
nonzeros = findall((output1 .> 0) .& (output2 .> 0))
sites = unique([nonzeros[i][2] for i in 1:length(nonzeros)])
betas = 1
