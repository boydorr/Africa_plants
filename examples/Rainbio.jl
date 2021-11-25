using JuliaDB
using JuliaDBMeta
using Unitful
using EcoSISTEM.ClimatePref

rainbio = loadtable("../../Data/RAINBIO/published_database/RAINBIO.csv")
ref = create_reference(0.75)
x = collect(JuliaDB.select(rainbio, :decimalLatitude))
y = collect(JuliaDB.select(rainbio, :decimalLongitude))
refval = extractvalues(y .* °, x .* °, ref)
rainbio = pushcol(rainbio, :refval, refval)
numspp = length(unique(JuliaDB.select(rainbio, :species)))

function startingArray(rainbio::JuliaDB.IndexedTable, numspecies::Int64, xmin = -25°, xmax = 50°, ymin = -35°, ymax = 40°)
    ref = create_reference(0.75)
    newref = ref.array[xmin.. xmax, ymin.. ymax]
    rainbio = filter(g->g.refval in newref, rainbio)
    refdict = Dict(zip(newref[1:end], 1:length(newref)))
    fillarray = Array{Int64, 2}(undef, numspecies, length(newref))
    ids = sort(unique(collect(JuliaDB.select(rainbio, :species))))
    dict = Dict(zip(ids, 1:length(ids)))
    grouped_tab = @groupby rainbio (:species, :refval) {count = length(:idrb)}
    sppnames = [dict[x] for x in collect(JuliaDB.select(grouped_tab, :species))]
    refs = [refdict[y] for y in collect(JuliaDB.select(grouped_tab, :refval))]
    counts = collect(JuliaDB.select(grouped_tab, :count))
    map(1:length(counts)) do i
        fillarray[sppnames[i], refs[i]] = counts[i]
    end
    return fillarray
end
africa = readfile("data/Africa.tif", -25°, 50°, -35°, 40°)[:, end:-1:1]
active =  Array{Bool, 2}(.!isnan.(africa))
start = startingArray(rainbio, numspp)
eco = Metacommunity(start)
diver = norm_sub_alpha(eco, 0)[:diversity]
diver = reshape(diver, 100, 100)
diver[isnan.(diver)] .= 0
diver[.!(active)] .= NaN
plotlyjs()
heatmap(-25:0.75:49.25, -35:0.75:39.25, transpose(diver), background_color = :lightblue, background_color_outside=:white, grid = false, color = cgrad(:algae_r, scale = :exp), aspect_ratio = 1, colorbar_title = "SR", size = (1200, 800), layout = (@layout [a b]), clim = (0, 2000))

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
diverA = catmean("data/Africa_full_", 10)
mydiver = diverA[1, :, :, 118, 1]
mydiver[isnan.(mydiver)] .= 0
mydiver[.!(active)] .= NaN
heatmap!(-25:0.75:49.25, -35:0.75:39.25, transpose(mydiver), background_color = :lightblue, background_color_outside=:white, grid = false, color = cgrad(:algae_r, scale = :exp), aspect_ratio = 1, colorbar_title = "SR", size = (1200, 800), layout = (@layout [a b]), clim = (0, 2000), subplot = 2)
Plots.pdf("plots/Rainbio.pdf")
