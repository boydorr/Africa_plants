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
    SOIL_LAYER_DEPTH

The depth of ERA5's top soil layer, 7 cm, over which a volumetric soil water fraction becomes a
depth of water, as EcoSISTEM's catalogue reads `swvl1`.
"""
const SOIL_LAYER_DEPTH = 7.0cm

"""
    species_pool(envelopes, widths; inflation = 0.1)

Return the species a run simulates, as a named tuple of `pool`, a `DataFrame` with one row per
species, and the names left out: `unmatched`, those with records but no envelope, and `nowidth`,
those whose records give no width on one of the two tolerance axes.

`pool` holds `species`, the temperature tolerance `tmin` and `tmin_sd` (K), the precipitation
tolerance `tp` and `tp_sd` (mm/day), and the envelope's mean exposures `ssr` (J m^-2 per day) and
`swvl1` (volumetric fraction), sorted by name. Each width has `inflation` times the pool's largest
width on that axis added, so that no species has a vanishing niche.

# Arguments

  - `envelopes`: the published envelopes, with `species`, `tmin`, `tp` (m/day), `ssr` and `swvl1`.
  - `widths`: the widths across each species' records, with `species`, `tmin_sd` (K) and `tp_sd`
    (m/day), `missing` where a species has too few records.
  - `inflation`: the share of the largest width added to every width.
"""
function species_pool(envelopes::DataFrame, widths::DataFrame;
                      inflation::Real = 0.1)
    known = Set(envelopes.species)
    unmatched = sort([s for s in widths.species if s ∉ known])
    joined = innerjoin(select(widths, :species, :tmin_sd, :tp_sd),
                       select(envelopes, :species, :tmin, :tp, :ssr, :swvl1),
                       on = :species)
    hasboth = .!ismissing.(joined.tmin_sd) .& .!ismissing.(joined.tp_sd)
    nowidth = sort(joined.species[.!hasboth])
    pool = joined[hasboth, :]
    tmin_sd = Float64.(pool.tmin_sd)
    tp_sd = Float64.(pool.tp_sd)
    pool.tmin_sd = tmin_sd .+ inflation * maximum(tmin_sd)
    pool.tp = pool.tp .* 1000
    pool.tp_sd = (tp_sd .+ inflation * maximum(tp_sd)) .* 1000
    sort!(pool, :species)
    return (pool = select(pool, :species, :tmin, :tmin_sd, :tp, :tp_sd, :ssr,
                          :swvl1),
            unmatched = unmatched, nowidth = nowidth)
end

"""
    envelopes_record(path)

Return the [`InputRecord`](@ref) of the published climate envelopes read from `path`, for a run's
provenance.
"""
function envelopes_record(path::AbstractString)
    return EcoSISTEM.InputRecord(role = :species,
                                 dataset = "Global plant climate envelopes",
                                 path = basename(path), doi = ENVELOPES_DOI,
                                 citation = "Harris, C., Brummitt, N., Cobbold, C. A. and Reeve, R. (2022). " *
                                            "NERC EDS Environmental Information Data Centre. " *
                                            "https://doi.org/$ENVELOPES_DOI")
end

"""
    africa_species(pool; density = 1.0 / m^2, survival = 0.2, seed = 1, provenance = EcoSISTEM.InputRecord[])

Build the [`SpeciesList`](@ref) of the species in `pool`, the `pool` of [`species_pool`](@ref),
named by their binomials and with no individuals.

Each species is matched to temperature and precipitation by a Gaussian tolerance of the pool's mean
and width, and demands three resources per individual: `1 / density` of ground, its mean net solar
radiation over that ground, and its mean soil water, as a depth over [`SOIL_LAYER_DEPTH`](@ref),
over that ground. Birth and death are 0.15 a year, longevity 1, and each species' mean dispersal
distance is drawn uniformly from 0.6 to 2.4 km.

# Arguments

  - `pool`: one row per species, as [`species_pool`](@ref) returns it.
  - `density`: plants per unit of available ground, which sets the area each individual occupies.
  - `survival`: how strongly a poor niche fit raises deaths and lowers births.
  - `seed`: seeds the dispersal distances and the species' own random streams.
  - `provenance`: the records of the data the pool was built from.
"""
function africa_species(pool::DataFrame; density::Unitful.Quantity = 1.0 / m^2,
                        survival::Real = 0.2, seed::Integer = 1,
                        provenance = EcoSISTEM.InputRecord[])
    n = nrow(pool)
    area = uconvert(m^2, 1 / density)
    dispersal = rand(Random.Xoshiro(seed), Uniform(0.6, 2.4), n) .* km
    perday = 1 / Unitful.d
    return build_species(n,
                         tolerance = (temperature = (pool.tmin .* K,
                                                     pool.tmin_sd .* K),
                                      rainfall = (pool.tp .* mm .* perday,
                                                  pool.tp_sd .* mm .* perday)),
                         toleranceaxis = (temperature = Temperature,
                                          rainfall = Precipitation),
                         demand = (space = fill(area, n),
                                   light = pool.ssr .* (J / m^2) .* perday .*
                                           area,
                                   water = pool.swvl1 .* SOIL_LAYER_DEPTH .*
                                           area),
                         demandaxis = (space = SurfaceArea,
                                       light = SolarRadiation,
                                       water = SoilWaterVolume),
                         dispersal = dispersal,
                         birth = 0.15 / year, death = 0.15 / year,
                         longevity = 1.0, survival = survival,
                         abundance = zeros(Int, n), names = pool.species,
                         seed = seed, provenance = provenance)
end
