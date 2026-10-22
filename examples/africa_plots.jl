# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Species richness across Africa from the forward runs: the first and last recorded years, the
# control's last year where control runs exist, and the departure from it, two panels to a row,
# averaged over every replicate seed in `AFRICA_SAVEDIR`, as the paper's richness figure is. The
# departure from the control is the part of the change the climate is responsible for; with no
# control runs to compare against, the change between the two years is drawn instead.
#
#     julia --project=examples examples/africa_plots.jl
#
# | variable         | default        | meaning                                         |
# |------------------|----------------|-------------------------------------------------|
# | `AFRICA_SAVEDIR` | `./Africa_run` | where `forward_seed*.jld2` are read and the figure written |
#
# Richness is normalised subcommunity alpha diversity at order 0, which for one cell is the number of
# species in it. An active cell with no individuals has no diversity to record and is drawn as zero;
# a cell outside the simulated land is left blank. Writes `richness.pdf` beside the runs.

using JLD2
using Plots
using Statistics

const SAVEDIR = get(ENV, "AFRICA_SAVEDIR", joinpath(pwd(), "Africa_run"))

# Richness runs from near black where nothing grows to pale green where everything does, as the
# published figure has it; the diverging panels keep their own scale.
const RICHNESS = cgrad(:algae, rev = true)

# Every run of `kind` in `SAVEDIR`, as the loaded files, refusing runs made on different grids.
function load_runs(kind::AbstractString)
    names = sort(filter(f -> startswith(f, kind * "_seed") &&
                             endswith(f, ".jld2"),
                        readdir(SAVEDIR)))
    runs = [load(joinpath(SAVEDIR, n)) for n in names]
    for r in runs
        haskey(r, "shape") ||
            error("$(kind) runs in $(SAVEDIR) do not record their grid; rerun them with the current " *
                  "`examples/africa_forward.jl`")
        (r["shape"] == runs[1]["shape"] && r["years"] == runs[1]["years"]) ||
            error("the $(kind) runs in $(SAVEDIR) are not all on one grid and span")
    end
    return runs
end

# Each cell's richness in every recorded year, averaged over `runs`: cells x years, zero in an empty
# active cell and `NaN` outside the simulated land.
function mean_richness(runs)
    q0 = findfirst(==(0.0), runs[1]["qs"])
    isnothing(q0) && error("the runs recorded no diversity at order 0")
    stack = [begin
                 d = r["diversity"][:, q0, :]
                 active = r["active"]
                 [active[c] ? (isnan(d[c, t]) ? 0.0 : d[c, t]) : NaN
                  for c in axes(d, 1), t in axes(d, 2)]
             end
             for r in runs]
    return reduce(+, stack) ./ length(stack)
end

# One year of `richness` as a `y x x` map on the runs' grid, rows in ascending `y` as `heatmap` expects.
function as_map(richness::AbstractMatrix, run, year)
    column = findfirst(==(year), run["years"])
    grid = reshape(richness[:, column], run["shape"])
    return run["y"][1] < run["y"][end] ? grid : grid[end:-1:1, :]
end

# A heatmap of `values` on the runs' grid in kilometres, with the continent's own aspect.
function panel(run, values; title, clims = :auto, color = RICHNESS)
    y = run["y"]
    ys = y[1] < y[end] ? y : reverse(y)
    return heatmap(run["x"] ./ 1000, ys ./ 1000, values, title = title,
                   titleloc = :left,
                   aspect_ratio = :equal, color = color, clims = clims,
                   background_color_inside = :lightblue, grid = false,
                   xlabel = "km", ylabel = "km",
                   xticks = -8000:4000:8000)
end

forward = load_runs("forward")
isempty(forward) && error("no forward runs in $(SAVEDIR)")
run = forward[1]
first_year, last_year = run["years"][1], run["years"][end]
richness = mean_richness(forward)
before = as_map(richness, run, first_year)
after = as_map(richness, run, last_year)
top = max(maximum(filter(!isnan, before)), maximum(filter(!isnan, after)))
change = after .- before
span = maximum(abs, filter(!isnan, change); init = 1.0)

panels = [panel(run, before, title = "A  $(first_year)", clims = (0, top)),
    panel(run, after, title = "B  $(last_year)", clims = (0, top))]

neutral = load_runs("neutral")
if isempty(neutral)
    push!(panels,
          panel(run, change, title = "C  change", clims = (-span, span),
                color = :balance))
else
    # The control's own richness is a map like the other two and shares their scale, so the three
    # can be read against each other; the climate's part is the difference in the fourth panel.
    control = as_map(mean_richness(neutral), neutral[1], last_year)
    departure = after .- control
    reach = maximum(abs, filter(!isnan, departure); init = 1.0)
    push!(panels,
          panel(run, control, title = "C  $(last_year) control",
                clims = (0, top)),
          panel(run, departure, title = "D  $(last_year) against the control",
                clims = (-reach, reach), color = :balance))
end

# Two panels to a row, so that the maps are as large as the page allows: three of them leave the
# fourth place empty, which is where the control's panel goes once control runs exist.
columns = min(2, length(panels))
rows = cld(length(panels), columns)
figure = plot(panels..., layout = (rows, columns),
              size = (560 * columns, 460 * rows),
              margin = 3Plots.mm, left_margin = 10Plots.mm)
out = joinpath(SAVEDIR, "richness.pdf")
savefig(figure, out)
println(length(forward), " forward run", length(forward) == 1 ? "" : "s",
        isempty(neutral) ? "" :
        " and $(length(neutral)) control run" *
        (length(neutral) == 1 ? "" : "s"),
        "; richness ", first_year, " to ", last_year, " drawn to ", out)
