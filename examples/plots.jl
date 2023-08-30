using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.ClimatePref
using JLD2
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
africa = readfile("data/Africa.tif", -25°, 50°, -35°, 40°)[:, end:-1:1]
active =  Array{Bool, 2}(.!isnan.(africa))

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
