using RCall
using DataFrames
using Plots
using StatsPlots
using LinearAlgebra
plotlyjs()
R"library(rdiversity)
qs <- seq(0, 20, 0.1)
pop <- matrix(c(1, 1), nrow = 2, ncol= 1)
div <-  qD(pop, qs)
gn <- meta_gamma(metacommunity(pop), qs)
gn$partition_name ='0% similar'
z <- matrix(c(1, 0.5, 0.5, 1), nrow =2, ncol=2)
z2 <- matrix(c(1, 0.75, 0.75, 1), nrow =2, ncol=2)
gz <- meta_gamma(metacommunity(pop, z), qs)
gz$partition_name = '50% similar'
gz2 <- meta_gamma(metacommunity(pop, z2), qs)
gz2$partition_name = '75% similar'
g <- rbind(gn, gz, gz2)"
@rget g

@df g plot(:q, :diversity, group = :partition_name, grid = false, xlab = "Order of diversity (q)", ylab = "Diversity", size = (800, 500), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, ylim = (1,2))
Plots.pdf("/Users/claireh/Documents/PhD/GIT/Thesis/introduction/figures/Diversityplot.pdf")

pop = zeros(20, 2)
pop[:, 1] .= 1
@rput pop
divvals = map(1:1:21) do i
    @rput i
    R"library(rdiversity)
    qs <- 1
    Z <- diag(20)
    Z[,] = 0.5
    diag(Z) = 1
    if (i > 1){
        pop[(i-1), 2] = 1
    }
    gn <- raw_meta_beta(metacommunity(pop), qs)
    gz <- raw_meta_beta(metacommunity(pop, Z), qs)
    gn$partition_name = (i-1)
    gz$partition_name = (i-1)
    "
    @rget gn
    @rget gz
    return gn, gz
end
div = vcat([divvals[x][1] for x in 1 :21]...)
div2 = vcat([divvals[x][2] for x in 1 :21]...)
@df div plot(:partition_name, :diversity, grid = false, xlab = "Number of shared species", ylab = "Distinctiveness", size = (800, 500), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, ylim = (0, 1), xlim= (0, 20), layout = (@layout [a b]))
@df div2 plot!(:partition_name, :diversity, color = 2, subplot = 1)
hline!([0.5], color = :black, linestyle = :dash, legend = false, subplot =1)

pop = zeros(20, 2)
pop[:, 1] .= 1
@rput pop
divvals = map(1:1:21) do i
    @rput i
    R"library(rdiversity)
    qs <- 1
    Z <- diag(20)
    Z[,] = 0.5
    diag(Z) = 1
    if (i > 1){
        pop[(i-1), 2] = 1
    }
    gn <- raw_meta_rho(metacommunity(pop), qs)
    gz <- raw_meta_rho(metacommunity(pop, Z), qs)
    gn$partition_name = (i-1)
    gz$partition_name = (i-1)
    "
    @rget gn
    @rget gz
    return gn, gz
end
div = vcat([divvals[x][1] for x in 1:21]...)
div2 = vcat([divvals[x][2] for x in 1:21]...)
@df div plot!(:partition_name, :diversity, grid = false, xlab = "Number of shared species", ylab = "Redundancy", size = (800, 500), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, ylim = (1, 2), xlim= (0, 20), subplot = 2, right_margin = 5.0*Plots.mm, color = 1)
@df div2 plot!(:partition_name, :diversity, color = 2, subplot = 2)
hline!([2.0], color = :black, linestyle = :dash, legend = false, subplot = 2)
Plots.pdf("/Users/claireh/Documents/PhD/GIT/Thesis/introduction/figures/Beta.pdf")

divvals = map(1:1:20) do i
    @rput i
    R"library(rdiversity)
    qs <- 1
    pop <- matrix(rep(1, i), nrow = i, ncol= 1)
    gn <- meta_gamma(metacommunity(pop), qs)
    gn$partition_name = i
    "
    @rget gn
    return gn
end
div = vcat(divvals...)
@df div plot(:partition_name, :diversity, grid = false, xlab = "Order of diversity (q)", ylab = "Diversity", size = (800, 500), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, ylim = (1, 20), xlim= (1, 20))
@df div plot!(:partition_name, log.(:diversity), ylim = (1, 20), xlim= (1, 20))

divvals = map(1:1:20) do i
    @rput i
    R"library(rdiversity)
    qs <- 1
    if (i == 0){
        pop <- matrix(1, nrow = 1, ncol= 1)
        gn <- meta_gamma(metacommunity(pop), qs)
        gn$diversity = 0
    }else{
        pop <- matrix(rep(i, i), nrow = i, ncol= 1)
        gn <- meta_gamma(metacommunity(pop), qs)
    }
    gn$partition_name = i
    "
    @rget gn
    return gn
end
div = vcat(divvals...)
@df div plot(:partition_name, :diversity, grid = false, xlab = "Number of species", ylab = "Value", size = (700, 500), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, ylim = (0, 20), xlim= (1, 20), aspect_ratio = 1, label = "Shannon diversity", layout = (@layout [a b]), link = :x)
#div[:diversity][1] = 1
@df div plot!(:partition_name, log.(:diversity), ylim = (0, 20), xlim= (1, 20), label = "Shannon entropy", subplot = 1, link = :x)

divvals = map(1:1:20) do i
    @rput i
    R"library(rdiversity)
    qs <- 2
    if (i == 0){
        pop <- matrix(1, nrow = 1, ncol= 1)
        gn <- meta_gamma(metacommunity(pop), qs)
        gn$diversity = 0
    }else{
        pop <- matrix(rep(i, i), nrow = i, ncol= 1)
        gn <- meta_gamma(metacommunity(pop), qs)
    }
    gn$partition_name = i
    "
    @rget gn
    return gn
end
div = vcat(divvals...)
@df div plot!(:partition_name, :diversity, grid = false, xlab = "Number of species", ylab = "", size = (700, 500), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, ylim = (0, 20), xlim= (1, 20), aspect_ratio = 1, label = "Simpson diversity", subplot = 2, color = 3, link = :x)
#div[:diversity][1] = 1
@df div plot!(:partition_name, 1 ./(:diversity), ylim = (0, 20), xlim= (1, 20), label = "Simpson index", subplot = 2, color = 4, link = :x)



Plots.pdf("/Users/claireh/Documents/PhD/GIT/Thesis/introduction/figures/Shannon.pdf")


using LaTeXStrings
pyplot()
ϵ = 0:0.01:1
plot(ϵ, ϵ.^-1, grid = false, legend = false, 
xlabel = "Species resource requirement,\$\\frac{\\epsilon}{\\bar\\epsilon}\$", 
ylabel = "Penalty for birth/boost for death, \$(\\frac{\\epsilon}{\\bar\\epsilon})^{-\\lambda}\$", 
guidefontsize = 18, tickfontsize = 16, size = (1600,600), linecolor = :black,
layout = (@layout [a b c]), top_margin = Plots.mm * 20.0,
annotations = [(0, 110, Plots.text("A", 24)),
(1.25, 110, Plots.text("B", 24)), (2.8, 110, Plots.text("C", 24))])

curr = 25
sds = 1:0.1:5.0
t = [(31, 0.0, Plots.text("Specialist", 16, rotation = 90)), 
(31, 0.9, Plots.text("Generalist", 16, rotation = 90))]
for sd in 1:length(sds)
    opt =  curr .+ (sds[sd] .* -5:0.05:5)
    match = (1.0)/sqrt(2 * pi * sds[sd]^2) .* exp.(-abs.(curr .- opt).^2 ./(2 * sds[sd]^2))
    plt = plot!(opt, match, subplot = 2, seriestype = :line, 
    legend = :none, xlims = (20, 30), color = :blues, 
    line_z = sds[sd], colorbar = :right, grid=nothing, 
    colorbar_title = "Niche width (°C)", 
    xlab = "Temperature (°C)", 
    ylab = "Match to environment, f", clim = (1, 5),
    guidefontsize = 18, tickfontsize = 16, ylim = (0, 1),
    annotations = t, right_margin = 20.0 * Plots.px)
end

curr = [25.5, 26.0, 26.5, 27.0, 27.5]
cols = cgrad([:gold, :orange, :darkorange, :darkred], [0.25, 0.5, 0.75, 1.0])
for i in eachindex(curr)
    vline!([curr[i]], linestyle = :dash, linewidth = 1, subplot = 2,
    color = cols[i * 1/5])
end
display(plt)

sds = collect(0.0:0.1:5.0)
opt =  25.0
for i in curr
    match = (1.0) ./ sqrt.(2 .* pi .* sds .^ 2) .* exp.(-abs.(i .- opt).^2 ./(2 .* sds.^2))
    plt = plot!(sds, match, seriestype = :line, legend = :topright, 
    grid=nothing, guidefontsize = 18, tickfontsize = 16,
    xlab = "Niche width (°C)", ylab = "Match to environment, f",
    legend_title = "Shift in temperature (°C)",
    subplot = 3, color = cgrad([:gold, :orange, :darkorange, :darkred], [0.25, 0.5, 0.75, 1.0]), 
    line_z = i - opt, label = i - opt, colorbar = :none)
end
display(plt)
Plots.pdf("plots/NicheWidth2.pdf")
