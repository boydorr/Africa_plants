# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Fetch the ERA5 monthly means the Africa runs use from the Climate Data Store: one file per
# variable per decade over the Africa box, through EcoSISTEM's `era5requests` and `fetchfiles`.
# A file already present is not asked for again, and EcoSISTEM writes a provenance record beside
# each file when it fetches it, or when it first finds it without one.
#
#     julia --project data/src/era5_download.jl [--dry-run]
#
# Needs a `~/.cdsapirc` holding the CDS URL and key, which CDSAPI.jl reads. Writes under
# `$AFRICA_DATADIR/era5/` (default `data/era5/`). A decade not yet complete is named by its years
# (`era5_t2m_2020-2025.nc`), so asking again once another year is complete fetches that block
# afresh; delete the old one.

module ERA5Download

using CDSAPI
using Dates
using EcoSISTEM

include("common.jl")

const CODES = ("t2m", "tp", "ssr", "swvl1")
const FIRST_YEAR = 1940

"""
    last_complete_year(today = today())

Return the last year every month of which is available at the CDS: the previous calendar year.
"""
last_complete_year(today::Date = Dates.today()) = Dates.year(today) - 1

"""
    cdsarea()

Return the Africa box as the CDS takes it, `(north, west, south, east)` in whole degrees.
"""
function cdsarea()
    return round.(Int,
                  (AFRICA_BOX.north, AFRICA_BOX.west, AFRICA_BOX.south,
                   AFRICA_BOX.east))
end

"""
    requests(dir; lastyear = last_complete_year())

Return, for each of `CODES` in order, the vector of `EcoSISTEM.CDSRequest`s that fetch that
variable from `FIRST_YEAR` to `lastyear` into `dir`, one per decade.
"""
function requests(dir::AbstractString;
                  lastyear::Integer = last_complete_year())
    return [EcoSISTEM.era5requests(code, FIRST_YEAR:lastyear, dir = dir,
                                   area = cdsarea())
            for code in CODES]
end

function main(argv)
    dir = joinpath(datadir(), "era5")
    mkpath(dir)
    specs = [SourceSpec(ERA, code, files = files)
             for (code, files) in zip(CODES, requests(dir))]
    if "--dry-run" in argv
        for entry in EcoSISTEM.fetchfiles(specs, dryrun = true)
            println(entry.present ? "  have   " : "  fetch  ",
                    basename(entry.entry.path))
        end
        return nothing
    end
    paths = EcoSISTEM.fetchfiles(specs)
    println(length(paths), " file(s) in ", dir)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ERA5Download.main(ARGS)
end
