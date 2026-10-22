# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The species pool: every species with seeding records in Africa that the published climate
# envelopes name, with the niche each is matched to, what it consumes and how it breeds and moves.

"""
    ENVELOPES_DOI

The DOI of the published plant climate envelopes the species' tolerance means and demands are
taken from.
"""
const ENVELOPES_DOI = "10.5285/ca339c86-3674-4030-b891-35326e71141e"

"""
    ENVELOPES_PAPER_DOI

The DOI of the paper the envelopes were published with, which describes how they were built.
"""
const ENVELOPES_PAPER_DOI = "10.1111/geb.13564"

"""
    PRIOR_RECORDS

The fewest records a species' width must rest on to inform the prior width of its relatives.
"""
const PRIOR_RECORDS = 30

"""
    PRIOR_MEMBERS

The fewest well-recorded relatives a taxonomic group must hold to give a prior width.
"""
const PRIOR_MEMBERS = 5

"""
    shrink_width(width, n, prior, strength)

Return the `width` across `n` records shrunk towards `prior`, which carries the weight of `strength`
further records: the two variances averaged by those weights. A width across few records is mostly
the prior, and one across many is mostly its own.
"""
function shrink_width(width::Real, n::Integer, prior::Real, strength::Real)
    return sqrt((n * width^2 + strength * prior^2) / (n + strength))
end

"""
    width_priors(widths, counts, taxonomy; records = PRIOR_RECORDS, members = PRIOR_MEMBERS)

Return the prior width of every species from its relatives, as a named tuple of the `prior` and the
`depth` of the level it came from. A species' prior is the median width of the well-recorded species
- those whose width rests on at least `records` records - in the first of its taxonomic groups to
hold at least `members` of them besides itself, its genus first; where none does, it is the median
over every well-recorded species, at a depth one past the last level.

  - `widths`: the width across each species' records, `missing` where it has too few.
  - `counts`: how many records each width rests on.
  - `taxonomy`: each species' groups from the narrowest outwards, such as `(genus, family, order)`,
    `missing` where one is not known.
"""
function width_priors(widths::AbstractVector, counts::AbstractVector{<:Integer},
                      taxonomy::AbstractVector{<:Tuple};
                      records::Integer = PRIOR_RECORDS,
                      members::Integer = PRIOR_MEMBERS)
    informed = [i
                for i in eachindex(widths)
                if !ismissing(widths[i]) && counts[i] >= records]
    # With nobody well recorded, every known width stands in, so that a small pool still has a prior.
    known = isempty(informed) ? findall(!ismissing, widths) : informed
    overall = median(Float64[widths[i] for i in known])
    depths = isempty(taxonomy) ? 0 : length(first(taxonomy))
    levels = [Dict{String, Vector{Int}}() for _ in 1:depths]
    for i in informed, (depth, name) in enumerate(taxonomy[i])
        ismissing(name) || push!(get!(levels[depth], name, Int[]), i)
    end
    prior = fill(overall, length(widths))
    depth = fill(depths + 1, length(widths))
    for i in eachindex(widths)
        for (d, name) in enumerate(taxonomy[i])
            ismissing(name) && continue
            others = Float64[widths[m]
                             for m in get(levels[d], name, Int[]) if m != i]
            length(others) >= members || continue
            prior[i] = median(others)
            depth[i] = d
            break
        end
    end
    return (prior = prior, depth = depth)
end

"""
    species_pool(envelopes, widths; strength = 12.0)

Return the species a run simulates, as a named tuple of `pool`, a `DataFrame` with one row per
species; `unmatched`, the names with records but no envelope, which are left out; and `nowidth`, the
names whose records give no width of their own on one of the two tolerance axes, which are kept on
their relatives' width.

`pool` holds `species`, the temperature tolerance `tmin` and `tmin_sd` (K), the precipitation
tolerance `tp` and `tp_sd` (mm/day), and the envelope's mean exposures `ssr` (J m^-2 per day) and
`swvl1` (volumetric fraction), sorted by name. `swvl1` is carried for comparison and is not used by
[`africa_species`](@ref), whose water demand is one rate for every species.

The width across a species' records says little when the records are few - it grows with their
number, and is zero for records from one place - so each width is shrunk towards that of the
species' well-recorded relatives by [`shrink_width`](@ref), the prior from [`width_priors`](@ref)
over genus, family and order. A species with too few records to have a width at all takes the prior
as it stands.

  - `envelopes`: the published envelopes, with `species`, `genus`, `family`, `order`, `tmin`, `tp`
    (m/day), `ssr` and `swvl1`.
  - `widths`: the widths across each species' records, with `species`, `tmin_sd` (K) and `tp_sd`
    (m/day), `missing` where a species has too few records, and the records each rests on, `tmin_n`
    and `tp_n`.
  - `strength`: how many records' weight the relatives' width carries against a species' own.
"""
function species_pool(envelopes::DataFrame, widths::DataFrame;
                      strength::Real = 12.0)
    known = Set(envelopes.species)
    unmatched = sort([s for s in widths.species if s ∉ known])
    joined = innerjoin(select(widths, :species, :tmin_n, :tmin_sd, :tp_n,
                              :tp_sd),
                       select(envelopes, :species, :genus, :family, :order,
                              :tmin, :tp, :ssr, :swvl1),
                       on = :species)
    taxonomy = collect(zip(joined.genus, joined.family, joined.order))
    own = (tmin_sd = copy(joined.tmin_sd), tp_sd = copy(joined.tp_sd))
    for (width, count) in ((:tmin_sd, :tmin_n), (:tp_sd, :tp_n))
        priors = width_priors(own[width], joined[!, count], taxonomy)
        joined[!, width] = Float64[ismissing(w) ? p :
                                   shrink_width(w, n, p, strength)
                                   for (w, n, p) in zip(own[width],
                                                        joined[!, count],
                                                        priors.prior)]
    end
    nowidth = sort(joined.species[ismissing.(own.tmin_sd) .| ismissing.(own.tp_sd)])
    pool = joined
    pool.tp = pool.tp .* 1000
    pool.tp_sd = pool.tp_sd .* 1000
    sort!(pool, :species)
    return (pool = select(pool, :species, :tmin, :tmin_sd, :tp, :tp_sd, :ssr,
                          :swvl1),
            unmatched = unmatched, nowidth = nowidth)
end

"""
    scaled_pool(pool, recorded)

Return the rows of `pool` a run at the scale `AFRICA_SCALE` names keeps, `recorded` holding the
species name of every seeding record. A small run keeps the most-recorded species, which fill a
coarse grid quickly. A medium run keeps species at even steps down the ranking by record count, so
that it holds the rarely recorded in the proportion the whole pool does. A large run keeps them all.
"""
function scaled_pool(pool::DataFrame, recorded::AbstractVector)
    (smallscale() || mediumscale()) || return pool
    counts = Dict{String, Int}()
    for s in recorded
        ismissing(s) || (counts[s] = get(counts, s, 0) + 1)
    end
    # By record count and then by name, so that the choice does not turn on the order of ties.
    ranked = sort(pool.species, by = s -> (-get(counts, s, 0), s))
    keep = smallscale() ? first(ranked, SCALED_SPECIES.small) :
           ranked[round.(Int,
                         range(1, length(ranked),
                               length = min(SCALED_SPECIES.medium,
                                            length(ranked))))]
    chosen = Set(keep)
    return filter(r -> r.species in chosen, pool)
end

"""
    envelopes_record(path)

Return the [`InputRecord`](@ref) of the published climate envelopes read from `path`, for a run's
provenance.
"""
function envelopes_record(path::AbstractString)
    return EcoSISTEM.InputRecord(role = :species,
                                 dataset = "Bioclimatic envelopes for global plant species",
                                 path = basename(path), doi = ENVELOPES_DOI,
                                 licence = "Open Government Licence v3",
                                 citation = "Harris, C., Brummitt, N., Cobbold, C. A. and Reeve, R. (2022). " *
                                            "Bioclimatic envelopes for global plant species. NERC EDS " *
                                            "Environmental Information Data Centre. " *
                                            "https://doi.org/$ENVELOPES_DOI. Described in Harris et al. " *
                                            "(2022), Strong phylogenetic signals in global plant " *
                                            "bioclimatic envelopes, Global Ecology and Biogeography 31, " *
                                            "2191-2203, https://doi.org/$ENVELOPES_PAPER_DOI")
end

"""
    africa_species(pool; density = 1.0 / m^2, water = nothing, survival = 0.2, seed = 1,
                   provenance = EcoSISTEM.InputRecord[])

Build the [`SpeciesList`](@ref) of the species in `pool`, the `pool` of [`species_pool`](@ref),
named by their binomials and with no individuals.

Each species is matched to temperature and precipitation by a Gaussian tolerance of the pool's mean
and width, and demands three resources per individual: `1 / density` of ground, its own mean net
solar radiation over that ground, and the rain falling on it at a rate that is the same for every
species. A common demand makes water limit how many plants a place holds without favouring one
species over another, which the precipitation tolerance already does. Birth and death are 0.15 a year, longevity 1, and each species' mean dispersal
distance is drawn uniformly from 0.6 to 2.4 km.

# Arguments

  - `pool`: one row per species, as [`species_pool`](@ref) returns it.
  - `density`: plants per unit of available ground, which sets the area each individual occupies.
  - `water`: the rainfall each individual demands over the ground it occupies, as a depth a day;
    `nothing` takes the median of the pool's precipitation optima, the rain a typical species is
    found under, so that a cell holds its full complement of plants only where it is at least that
    wet.
  - `survival`: how strongly a poor niche fit raises deaths and lowers births.
  - `seed`: seeds the dispersal distances and the species' own random streams.
  - `provenance`: the records of the data the pool was built from.
"""
function africa_species(pool::DataFrame; density::Unitful.Quantity = 1.0 / m^2,
                        water::Union{Nothing, Unitful.Quantity} = nothing,
                        survival::Real = 0.2, seed::Integer = 1,
                        provenance = EcoSISTEM.InputRecord[])
    n = nrow(pool)
    area = uconvert(m^2, 1 / density)
    dispersal = rand(Random.Xoshiro(seed), Uniform(0.6, 2.4), n) .* km
    perday = 1 / Unitful.d
    rain = isnothing(water) ? median(pool.tp) * mm * perday : water
    return build_species(n,
                         tolerance = (temperature = (pool.tmin .* K,
                                                     pool.tmin_sd .* K),
                                      rainfall = (pool.tp .* mm .* perday,
                                                  pool.tp_sd .* mm .* perday)),
                         toleranceaxis = (temperature = AnnualMinimumTemperature,
                                          rainfall = AnnualMeanPrecipitation),
                         demand = (space = fill(area, n),
                                   light = pool.ssr .* (J / m^2) .* perday .*
                                           area,
                                   water = fill(uconvert(Unitful.L / Unitful.d,
                                                         rain * area), n)),
                         demandaxis = (space = SurfaceArea,
                                       light = SolarRadiation,
                                       water = Precipitation),
                         dispersal = dispersal,
                         birth = 0.15 / year, death = 0.15 / year,
                         longevity = 1.0, survival = survival,
                         abundance = zeros(Int, n), names = pool.species,
                         seed = seed, provenance = provenance)
end
