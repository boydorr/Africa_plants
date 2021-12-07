using JuliaDB
using JuliaDBMeta
using EcoSISTEM.ClimatePref
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM
using EcoSISTEM.Units
using AxisArrays
using JLD2
using Printf
import EcoSISTEM: SavedLandscape, HabitatUpdate

function emptypopulate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AB, rel::R) where {AB <: EcoSISTEM.AbstractAbiotic, R <: EcoSISTEM.AbstractTraitRelationship}
@warn "Ecosystem not populated!"
end


function reverseBurnin!(eco::Ecosystem, duration::Unitful.Time, timestep::Unitful.Time, cacheInterval::Unitful.Time, cacheFolder::String, gbif::JuliaDB.DIndexedTable, effort::Array{Float64,2}, fillarray::Array{Int64,2}, xmin = -25°, xmax = 50°, ymin = -35°, ymax = 40°)
    reverseTime!(eco)
    gbif = @transform gbif {intdate = round(Int64, uconvert.(month, :date .- 1901year) ./ month)}
    gbif = reindex(gbif, :intdate)
    ref = create_reference(0.75)
    newref = ref.array[xmin.. xmax, ymin.. ymax]
    refdict = Dict(zip(newref[1:end], 1:length(newref)))
    times = length(0s:timestep:duration)
    effort = hcat(effort[1:end]...)
    for i in 1:times
        if i in collect(select(gbif, :intdate))
            eco.abundances.matrix .+= round.(Int64, db_to_array(gbif, size(eco.abundances.matrix, 1), i, fillarray, refdict) .* effort)
        end
        update!(eco, timestep)
        if mod(i*timestep, cacheInterval) == 0.0year
            JLD.save(joinpath(cacheFolder, (@sprintf "%03d.jld" times/12)), "gl", SavedLandscape(eco.abundances))
        end
    end
    forwardTime!(eco)
end

function reverseTime!(eco::Ecosystem)
    newhabchange1 = HabitatUpdate(reverseEraChange, eco.abenv.habitat.h1.change.changefun)
    newhabchange2 = HabitatUpdate(reverseEraChange, eco.abenv.habitat.h2.change.changefun)
    eco.abenv.habitat.h1.change = newhabchange1
    eco.abenv.habitat.h2.change = newhabchange2
    eco.abenv.habitat.h1.time = size(eco.abenv.habitat.h1.matrix, 3)
    eco.abenv.habitat.h2.time = size(eco.abenv.habitat.h2.matrix, 3)
end
function forwardTime!(eco::Ecosystem)
    newhabchange1 = HabitatUpdate(eraChange, eco.abenv.habitat.h1.change.changefun)
    newhabchange2 = HabitatUpdate(eraChange, eco.abenv.habitat.h2.change.changefun)
    eco.abenv.habitat.h1.change = newhabchange1
    eco.abenv.habitat.h2.change = newhabchange2
    eco.abenv.habitat.h1.time = 1
    eco.abenv.habitat.h2.time = 1
end

function reverseEraChange(eco::Ecosystem, hab::ContinuousTimeHab, timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    hab.time -= round(Int64, monthstep/month)
    if hab.time < 1
        hab.time = size(hab.matrix, 3)
        @warn "More timesteps than available, have repeated"
    end
end

function db_to_array(gbif::JuliaDB.DIndexedTable, numspecies::Int64, time::Int64, fillarray::Array{Int64, 2}, refdict::Dict)
    gbif = filter(g-> (g.intdate == time) & (g.refval in keys(refdict)), gbif)
    ids = sort(unique(collect(select(gbif, :SppID))))
    dict = Dict(zip(ids, 1:length(ids)))
    grouped_tab = @groupby gbif (:SppID, :refval) {count = length(:UID)}
    sppnames = [dict[x] for x in collect(select(grouped_tab, :SppID))]
    refs = [refdict[y] for y in collect(select(grouped_tab, :refval))]
    counts = collect(select(grouped_tab, :count))
    map(1:length(counts)) do i
        fillarray[sppnames[i], refs[i]] = counts[i]
    end
    return fillarray
end


function startingArray(gbif::JuliaDB.DIndexedTable, numspecies::Int64, xmin = -25°, xmax = 50°, ymin = -35°, ymax = 40°)
    ref = create_reference(0.75)
    newref = ref.array[xmin.. xmax, ymin.. ymax]
    gbif_fil = filter(g->g.refval in newref, gbif)
    refdict = Dict(zip(newref[1:end], 1:length(newref)))
    fillarray = Array{Int64, 2}(undef, numspecies, length(newref))
    ids = sort(unique(collect(select(gbif, :SppID))))
    dict = Dict(zip(ids, 1:length(ids)))
    grouped_tab = @groupby gbif_fil (:SppID, :refval) {count = length(:UID)}
    sppnames = [dict[x] for x in collect(select(grouped_tab, :SppID))]
    refs = [refdict[y] for y in collect(select(grouped_tab, :refval))]
    counts = collect(select(grouped_tab, :count))
    map(1:length(counts)) do i
        fillarray[sppnames[i], refs[i]] = counts[i] .* 1e3
    end
    return fillarray
end

function startingArray(gbif::JuliaDB.DIndexedTable, numspecies::Int64, kernel::Bool, xmin = -25°, xmax = 50°, ymin = -35°, ymax = 40°)
    ref = create_reference(0.75)
    newref = ref.array[xmin.. xmax, ymin.. ymax]
    gbif_fil = filter(g->g.refval in newref, gbif)
    refdict = Dict(zip(newref[1:end], 1:length(newref)))
    fillarray = Array{Int64, 2}(undef, numspecies, length(newref))
    ids = sort(unique(collect(select(gbif, :SppID))))
    dict = Dict(zip(ids, 1:length(ids)))
    grouped_tab = @groupby gbif_fil (:SppID, :refval) {count = length(:UID)}
    sppnames = [dict[x] for x in collect(select(grouped_tab, :SppID))]
    refs = [refdict[y] for y in collect(select(grouped_tab, :refval))]
    counts = collect(select(grouped_tab, :count))
    neighbours = [get_neighbours(zeros(100, 100), convert_coords(refs[i], 100)[1],convert_coords(refs[i], 100)[2], 8) for i in eachindex(refs)]
    neighbours = [convert_coords.(neighbours[x][:, 1], neighbours[x][:, 2], 100) for x in eachindex(neighbours)]
    map(1:length(counts)) do i
        fillarray[sppnames[i], refs[i]] = counts[i]
        fillarray[sppnames[i], neighbours[i]] .= Int64(1e6)
    end
    return fillarray
end

# using GLM
# using DataFrames
# function linmod(x)
#     df = DataFrame(X = x, Y = 1:length(x))
#     mod = GLM.lm(@formula(X ~ Y), df)
#     return coef(mod)[2]
# end
#
#
# function slopediv(mat::Array{Float64, 5})
#     return mapslices(linmod, mat, dims = 4)[:, :, :, 1, :]
# end
# function slopediv(mat::Array{Float64, 6})
#     return mapslices(linmod, mat, dims = 4)[:, :, :, 1, :, :]
# end
