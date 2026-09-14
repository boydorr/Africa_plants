# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Fetch the 20CRv3 monthly means the extension before 1940 uses from NOAA PSL: one global netCDF
# file per variable, 1806 to 2015, over plain HTTPS with no credentials. The URLs come from
# EcoSISTEM's `TwentyCR` catalogue and the transfer from its `fetchfiles`, which resumes an
# interrupted file and writes a provenance record beside each one.
#
#     julia --project data/src/twentycr_download.jl [--dry-run]
#
# Writes under `$AFRICA_DATADIR/twentycr/SI-MO/`, keeping PSL's file names. PSL serves the same
# files under its `SI` directories too. The files are global and are not cut here.

module TwentyCRDownload

using EcoSISTEM

include("common.jl")

const CODES = ("air", "prate", "dswrf", "uswrf", "soilw")

"""
    assets(dir)

Return an `EcoSISTEM.CachedAsset` for each of `CODES` in order: the layer's download URL from the
`TwentyCR` catalogue, kept at `dir/SI-MO/` under the file name PSL gives it.
"""
function assets(dir::AbstractString)
    return [begin
                url = EcoSISTEM.layerinfo(TwentyCR, code).file
                EcoSISTEM.CachedAsset(TwentyCR, url,
                                      path = joinpath(dir, "SI-MO",
                                                      basename(url)))
            end
            for code in CODES]
end

function main(argv)
    dir = joinpath(datadir(), "twentycr")
    mkpath(joinpath(dir, "SI-MO"))
    specs = [SourceSpec(TwentyCR, code, file = asset)
             for (code, asset) in zip(CODES, assets(dir))]
    if "--dry-run" in argv
        for entry in EcoSISTEM.fetchfiles(specs, dryrun = true)
            println(entry.present ? "  have   " : "  fetch  ",
                    basename(entry.entry.path),
                    isnothing(entry.bytes) ? "" :
                    string(" (", round(entry.bytes / 2^20, digits = 1),
                           " MiB)"))
        end
        return nothing
    end
    paths = EcoSISTEM.fetchfiles(specs)
    println(length(paths), " file(s) in ", joinpath(dir, "SI-MO"))
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    TwentyCRDownload.main(ARGS)
end
