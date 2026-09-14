# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Fetch the country centroids that `gbif_clean.jl` excludes records near: the points a record's
# coordinates were set to when all that was known of its locality was the country's name.
#
#     julia --project data/src/centroids.jl [--force]
#
# Reads GBIF's Catalogue of Centroids at a pinned commit (centroids gathered from
# CoordinateCleaner, GeoNames, the Getty Thesaurus of Geographic Names, GEOLocate and others)
# and keeps its country-level (`PCLI`) entries. Province centroids and capitals are listed there
# too and are not kept: records near them are mostly at real localities. Writes `centroids.csv`
# with `iso2, name, source, latitude, longitude`. The catalogue states no licence, so the table
# is regenerated rather than committed.

module Centroids

using CSV
using DataFrames
using Downloads

include("common.jl")

const COMMIT = "0959e14d0e7a90de28ae2bc23cf570b791e112ef"
const URL = "https://raw.githubusercontent.com/gbif/catalogue-of-centroids/$COMMIT/centroids.tsv"
const TYPES = ["PCLI"]

"""
    country_centroids(tsv; types = TYPES)

Return the centroids in the Catalogue of Centroids file `tsv` whose `type` is one of `types`,
as a `DataFrame` with columns `iso2, name, source, latitude, longitude`; rows without
coordinates are dropped. Namibia's ISO code `NA` is kept as a code, not read as missing.
"""
function country_centroids(tsv::AbstractString; types = TYPES)
    raw = CSV.read(tsv, DataFrame, delim = '\t', missingstring = "",
                   types = Dict(:iso2 => String))
    rows = filter(r -> !ismissing(r.type) && r.type in types &&
                       !ismissing(r.decimal_latitude) &&
                       !ismissing(r.decimal_longitude), raw)
    return DataFrame(iso2 = coalesce.(rows.iso2, ""),
                     name = coalesce.(rows.gbif_name, ""),
                     source = coalesce.(rows.source, ""),
                     latitude = Float64.(rows.decimal_latitude),
                     longitude = Float64.(rows.decimal_longitude))
end

function main(argv)
    out = joinpath(datadir(), "centroids.csv")
    if isfile(out) && !forced(argv)
        println("Already have ", out,
                "; nothing to do (use --force to rebuild).")
        return nothing
    end
    tmp = out * ".part"
    println("Fetching ", URL, " ...")
    Downloads.download(URL, tmp)
    table = country_centroids(tmp)
    rm(tmp)
    CSV.write(out, table)
    write_provenance(out, source = "GBIF Catalogue of Centroids, centroids.tsv",
                     url = URL, commit = COMMIT, types = TYPES,
                     centroids = nrow(table),
                     citation = "GBIF Catalogue of Centroids, https://github.com/gbif/catalogue-of-centroids (no licence stated)")
    println(nrow(table), " centroids -> ", out)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Centroids.main(ARGS)
end
