# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Fetch the botanic gardens that `gbif_clean.jl` excludes records near: a record at a garden may be
# a cultivated plant, or a specimen given the coordinates of the garden that holds it, and says
# nothing about where the species grows wild.
#
#     julia --project data/src/gardens.jl [--force]
#
# Reads CoordinateCleaner's `institutions` gazetteer at a pinned commit (biodiversity institutions
# compiled from GBIF, Wikipedia, GeoNames, the Global Registry of Biodiversity Repositories, Index
# Herbariorum and Botanic Gardens Conservation International) and keeps every botanic garden in it:
# the rows typed as one, and the rows of any other type whose name says it is one. Writes
# `gardens.csv` with `name, type, country, source, latitude, longitude`. CoordinateCleaner is GPL-3,
# so the table is regenerated rather than committed.

module Gardens

# CodecBzip2 must be loaded before RData reads the file, which R saved bzip2-compressed.
using CodecBzip2
using CSV
using DataFrames
using Downloads
using RData

include("common.jl")

const COMMIT = "24706b5d9ac928fd3331a0cff9572caa7f575885"
const URL = "https://raw.githubusercontent.com/ropensci/CoordinateCleaner/$COMMIT/data/institutions.rda"

# A name that says the institution is a botanic garden or an arboretum, in the languages the
# gazetteer's names are written in. Deliberately narrower than "botan" or "garden", which would take
# in botany departments, herbaria and botanical societies.
const GARDEN_NAME = r"botanic(al)? garden|botanische[rn]? g(a|ä)rten|jard[ií]n bot|jardim bot|jardin des plantes|hortus bot|orto botanico|giardino botanico|ogr(ó|o)d botaniczny|arboret"i

# Whether a coordinate is present and a number.
_known(x) = !ismissing(x) && isfinite(x)

"""
    botanic_gardens(institutions)

Return the botanic gardens among the rows of CoordinateCleaner's `institutions` table that have
coordinates: those whose `type` is `Botanic_garden`, and those of any other type whose `name`
matches `GARDEN_NAME`. The result has columns `name, type, country, source, latitude, longitude`.
"""
function botanic_gardens(institutions::DataFrame)
    typed = coalesce.(institutions.type .== "Botanic_garden", false)
    named = [!ismissing(n) && occursin(GARDEN_NAME, n)
             for n in institutions.name]
    located = _known.(institutions.decimalLatitude) .&
              _known.(institutions.decimalLongitude)
    rows = institutions[(typed .| named) .& located, :]
    return DataFrame(name = coalesce.(rows.name, ""),
                     type = coalesce.(rows.type, ""),
                     country = coalesce.(rows.country, ""),
                     source = coalesce.(rows.source, ""),
                     latitude = Float64.(rows.decimalLatitude),
                     longitude = Float64.(rows.decimalLongitude))
end

function main(argv)
    out = joinpath(datadir(), "gardens.csv")
    if isfile(out) && !forced(argv)
        println("Already have ", out,
                "; nothing to do (use --force to rebuild).")
        return nothing
    end
    # RData picks its reader by the file's extension, so the download keeps its own name.
    tmp = joinpath(mktempdir(), "institutions.rda")
    println("Fetching ", URL, " ...")
    Downloads.download(URL, tmp)
    institutions = RData.load(tmp)["institutions"]
    table = botanic_gardens(institutions)
    CSV.write(out, table)
    typed = count(==("Botanic_garden"), table.type)
    write_provenance(out, source = "CoordinateCleaner institutions gazetteer",
                     url = URL, commit = COMMIT, gardens = nrow(table),
                     typed_botanic_garden = typed,
                     named_as_garden = nrow(table) - typed,
                     licence = "GPL-3 (CoordinateCleaner)",
                     citation = "Zizka et al. (2019) CoordinateCleaner: standardized cleaning of occurrence records from biological collection databases. Methods in Ecology and Evolution 10:744-751, doi:10.1111/2041-210X.13152")
    println(nrow(table), " botanic gardens (", typed, " typed as gardens, ",
            nrow(table) - typed, " by name) -> ", out)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Gardens.main(ARGS)
end
