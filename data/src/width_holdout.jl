# SPDX-License-Identifier: LGPL-3.0-or-later
#
# How well can a species' niche width be recovered from a few records? Species with many records,
# whose width across all of them is taken as known, are cut down to a few records each and their
# width estimated again by each candidate rule, and the estimates scored against the known width.
#
#     julia --project -t auto data/src/width_holdout.jl
#
# A rule shrinks the width across the few records towards a prior width with the weight of `k`
# further records, by the package's own `shrink_width`. The prior is either
# the `pool`'s, the median width of every well-recorded species, or the `taxon`'s, the median width
# of the well-recorded species of the same genus, or of the family or order where the genus holds too
# few, the species under test always left out, from the package's `width_priors`. The rule the models have used so far, a tenth of the
# pool's largest width added to every width, is scored beside them as `inflation`.
#
# The few records are drawn three ways, since the records of a rarely collected species come from a
# smaller area than a random few of a common one's, though not from one spot: `random` from all the
# species' records; `nearest`, those nearest one drawn at random; and `matched`, from the records
# within a distance of one drawn at random, the distance drawn from how far apart the records of the
# species that really have that few lie. Widths are those of `species_widths.jl`: the spread across records of the minimum
# monthly temperature and of the mean precipitation over the twelve months to each record, with the
# same records left out. Writes `width_holdout.csv`, one row per axis, sampling, record count, prior
# and `k`.

module WidthHoldout

using CSV
using DataFrames
using Random
using Statistics

using Africa: PRIOR_MEMBERS, PRIOR_RECORDS, shrink_width, width_priors

include("species_widths.jl")
using .SpeciesWidths: record_windows
using .SpeciesWidths.SpeciesTraits: load_stack, datadir, write_sidecar

# A species' width is taken as known from this many records.
const KNOWN_RECORDS = 100
const SAMPLE_SIZES = (2, 3, 5, 10)
const STRENGTHS = (0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 50.0,
                   100.0,
                   1000.0)
const ENVELOPES = get(ENV, "AFRICA_ENVELOPES",
                      expanduser("~/git/plant_trait/data/global_plant_climate_envelopes.csv"))

# The distance over the ground between two records, in km, flat at the scale of a species' records.
function _apart(lat1, lon1, lat2, lon2)
    return 111.2 * sqrt((lat1 - lat2)^2 + (cosd(lat1) * (lon1 - lon2))^2)
end

# The greatest distance between any two of the records at `lat` and `lon`, in km.
function _extent(lat::AbstractVector{Float64}, lon::AbstractVector{Float64})
    return maximum(_apart(lat[a], lon[a], lat[b], lon[b])
                   for a in eachindex(lat), b in eachindex(lat))
end

# The `n` records of a species to estimate its width from. `:random` draws them from all its
# records; `:nearest` takes those nearest one drawn at random; `:matched` draws them from the records
# within `reach` km of one drawn at random, the nearest standing in where too few are within it.
function _sample(rng, lat::Vector{Float64}, lon::Vector{Float64}, n::Int,
                 scheme::Symbol, reach::Float64)
    scheme == :random && return randperm(rng, length(lat))[1:n]
    seed = rand(rng, eachindex(lat))
    away = [_apart(lat[seed], lon[seed], lat[j], lon[j])
            for j in eachindex(lat)]
    within = findall(<=(reach), away)
    (scheme == :nearest || length(within) < n) &&
        return partialsortperm(away, 1:n)
    return within[randperm(rng, length(within))[1:n]]
end

# How far estimates lie from the known widths, as the natural log of their ratio: its root mean
# square and median size over the estimates that are not zero, and how many were.
function _score(estimates::Vector{Float64}, known::Vector{Float64})
    ratios = [log(e / t) for (e, t) in zip(estimates, known) if e > 0]
    return (rmse = sqrt(mean(abs2, ratios)), median_abs = median(abs.(ratios)),
            bias = exp(median(ratios)), zero = count(==(0), estimates))
end

"""
    holdout(values, lat, lon, species, taxonomy; rng)

Score every rule on one axis. `values`, `lat` and `lon` hold the axis's trusted records and
`species` the name of each; `taxonomy` maps a name to its `(genus, family, order)`. Return the
table of scores and how many species were tested.
"""
function holdout(values::Vector{Float64}, lat::Vector{Float64},
                 lon::Vector{Float64}, species::Vector{String},
                 taxonomy::Dict{String, NTuple{3, String}};
                 rng = Xoshiro(1))
    rows = Dict{String, Vector{Int}}()
    for (r, s) in enumerate(species)
        haskey(taxonomy, s) && push!(get!(rows, s, Int[]), r)
    end
    # The well-recorded species, which give the priors, in a fixed order.
    informed = sort([s for (s, rs) in rows if length(rs) >= PRIOR_RECORDS])
    widths = [std(values[rows[s]]) for s in informed]
    names = [taxonomy[s] for s in informed]
    inflation = 0.1 * maximum(std(values[rs])
                        for rs in Base.values(rows)
                        if length(rs) >= 2)
    tested = [i
              for (i, s) in enumerate(informed)
              if length(rows[s]) >= KNOWN_RECORDS && widths[i] > 0]
    known = widths[tested]
    # One species in tens of thousands does not move the pool's median, so it is not left out.
    pool = fill(median(widths), length(tested))
    related = width_priors(widths, [length(rows[s]) for s in informed], names)
    taxon = related.prior[tested]
    depths = related.depth[tested]
    println("  ", length(informed), " well-recorded species give priors; ",
            length(tested), " tested; taxon prior from genus ",
            count(==(1), depths), ", family ", count(==(2), depths), ", order ",
            count(==(3), depths), ", pool ", count(==(4), depths))
    flush(stdout)
    scores = DataFrame(sampling = String[], n = Int[], prior = String[],
                       k = Float64[], rmse = Float64[], median_abs = Float64[],
                       bias = Float64[], zero = Int[])
    for scheme in (:random, :matched, :nearest), n in SAMPLE_SIZES
        # How far apart the records of the species that really have `n` of them lie, which the
        # matched samples are drawn to resemble.
        real = [_extent(lat[rs], lon[rs])
                for rs in Base.values(rows)
                if length(rs) == n]
        picks = [begin
                     rs = rows[informed[i]]
                     rs[_sample(rng, lat[rs], lon[rs], n, scheme,
                                rand(rng, real))]
                 end
                 for i in tested]
        few = [std(values[pick]) for pick in picks]
        println("  ", scheme, " n = ", n, ": median extent ",
                round(Int,
                      median(_extent(lat[pick], lon[pick]) for pick in picks)),
                " km against ", round(Int, median(real)), " km for the ",
                length(real), " species with ", n, " records")
        flush(stdout)
        sampling = String(scheme)
        for k in STRENGTHS,
            (label, priors) in (("pool", pool),
            ("taxon", taxon))
            est = [shrink_width(s, n, p, k) for (s, p) in zip(few, priors)]
            push!(scores, (sampling, n, label, k, _score(est, known)...))
        end
        push!(scores,
              (sampling, n, "inflation", NaN,
               _score(few .+ inflation, known)...))
    end
    return (scores = scores, tested = length(tested))
end

function main()
    records = CSV.read(joinpath(datadir(), "records_africa.csv"), DataFrame,
                       select = [:species, :decimalLatitude, :decimalLongitude,
                           :year, :month, :on_lattice, :near_garden,
                           :garden_named])
    filter!(r -> !ismissing(r.species), records)
    println("Read ", nrow(records), " named records")
    flush(stdout)
    envelopes = CSV.read(ENVELOPES, DataFrame,
                         select = [:species, :genus, :family, :order])
    taxonomy = Dict(String(r.species) => (String(r.genus), String(r.family),
                     String(r.order))
                    for r in eachrow(dropmissing(envelopes)))
    era5 = joinpath(datadir(), "era5")
    windows = record_windows(records, load_stack(era5, "t2m"),
                             load_stack(era5, "tp"))
    trusted = (temperature = .!records.on_lattice,
               rainfall = .!(records.on_lattice .| records.near_garden .|
               records.garden_named))
    # Precipitation in mm/day, as the models hold it, so that its widths read as theirs do.
    values = (temperature = windows.tmin, rainfall = windows.tp .* 1000)
    tables = DataFrame[]
    for axis in (:temperature, :rainfall)
        keep = trusted[axis] .& .!isnan.(values[axis])
        println(axis, ": ", count(keep), " trusted records")
        result = holdout(values[axis][keep],
                         Float64.(records.decimalLatitude[keep]),
                         Float64.(records.decimalLongitude[keep]),
                         String.(records.species[keep]), taxonomy)
        insertcols!(result.scores, 1, :axis => String(axis))
        push!(tables, result.scores)
    end
    out = joinpath(datadir(), "width_holdout.csv")
    scores = vcat(tables...)
    CSV.write(out, scores)
    write_sidecar(out,
                  source = "records_africa.csv joined to ERA5 monthly means, with the taxonomy of the published envelopes",
                  known_from_records = KNOWN_RECORDS,
                  prior_from_records = PRIOR_RECORDS,
                  prior_members = PRIOR_MEMBERS,
                  scores = "natural log of estimate over known width: rmse, median absolute value, exp(median) as bias; zero counts estimates of zero, left out of the rest")
    println(nrow(scores), " scores -> ", out)
    return scores
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    WidthHoldout.main()
end
