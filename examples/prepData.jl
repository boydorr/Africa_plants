using JuliaDB
using JuliaDBMeta
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using JLD
using Statistics
using StatsBase
using OnlineStats

function adjust(x, adj, min, max)
    edges = range(min, stop = max, length = 1000)
    h = Hist(edges)
    fit!(h, x)
    counts = h.counts .* adj
    step = edges[2]-edges[1]
    edges = collect(edges) .+ step/2
    return mean(edges[1:(end-1)], weights(counts))
end


gbif = JuliaDB.load("gbif/full_data/Small_GBIF_new")
continents = JuliaDB.load("../sdb/PHYLO/Continents")
continents = distribute(continents, 1)

gbif_new = join(gbif, continents, how = :left, rkey = :refval, lkey = :refval)
gbif_new = filter(t-> !ismissing(t.continent), gbif_new)
gbif_new = filter(t-> t.continent == 1.0, gbif_new)
JuliaDB.save(gbif_new, "Chapter5/data/Full_GBIF_africa_new")

traits = JuliaDB.load("../sdb/PHYLO/CERA_JOIN_SIMPLE")
continents = JuliaDB.load("../sdb/PHYLO/Continents")
continents = distribute(continents, 1)
traits_new = join(traits, continents, how = :left, rkey = :refval, lkey = :refval)
traits_new = filter(t-> !ismissing(t.continent), traits_new)
traits_new = filter(t-> t.continent == 1.0, traits_new)
evi_counts = JLD.load("../sdb/PHYLO/Total_evi_counts_1.jld", "total")
gbif_counts = JLD.load("../sdb/PHYLO/Total_gbif_counts_1.jld", "total")
adjustment = gbif_counts ./ evi_counts
adjustment[isnan.(adjustment)] .= 1
adjustment[isinf.(adjustment)] .= 1
mins = [197.0K, 197.0K, 197.0K, 0K, 197.0K, 197.0K, 197.0K, 197.0K, 0.0m^3, 0.0m^3, 0.0m^3, 0.0m^3, 0.0J/m^2, 0.0m]
maxs = [320.0K, 320.0K, 320.0K, 80K, 320.0K, 320.0K, 320.0K, 320.0K, 1.0m^3, 1.0m^3, 1.0m^3, 1.0m^3, 3.0e7J/m^2, 0.1m]

phylo_traits_adj = @groupby traits_new :SppID {tmin_mean = adjust(uconvert.(K, :tmin), adjustment[:, 1], mins[1], maxs[1]),
    tmin_sd = std(uconvert.(K, :tmin)), 
    swvl1 = adjust(:swvl1mean, adjustment[:, 9], mins[9], maxs[9]),
    ssr =  adjust(:ssrmean, adjustment[:, 13], mins[13], maxs[13]), 
    tp_mean =  adjust(:tpmean, adjustment[:, 14], mins[14], maxs[14]),
    tp_sd = std(:tpmean)}
JuliaDB.save(phylo_traits_adj, "Chapter5/data/Africa_traits_new")


gbif = JuliaDB.load("data/Full_GBIF_africa_new")
traits = JuliaDB.load("data/Africa_traits_new")
traits = filter(tr -> (tr.swvl1 > 0.0m^3) & (tr.ssr > 0.0J/m^2) & (tr.tp_sd > 0.0mm) & (tr.tmin_sd > 0.0K), traits)
gbif = filter(gb -> gb.SppID in collect(select(traits, :SppID)), gbif)
JuliaDB.save(gbif, "data/GBIF_africa_fil")
JuliaDB.save(traits, "data/Africa_traits_fil")