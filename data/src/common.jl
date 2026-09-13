# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Shared by the data scripts: where the data lives, the Africa box, and provenance sidecars.
# Included by each script; not part of the `Africa` module, because none of this is needed by a
# run once the data exist.

using Dates
using TOML

"""
    datadir()

Return the directory the data scripts read and write, `AFRICA_DATADIR` if set and otherwise
`data/` at the repository root, created if missing.
"""
function datadir()
    dir = abspath(get(ENV, "AFRICA_DATADIR", joinpath(@__DIR__, "..")))
    mkpath(dir)
    return dir
end

"""
    AFRICA_BOX

The geographic box every download is cut to, in degrees: the continent with Madagascar and a
margin of about a degree. The same numbers as `era5_request.toml`'s `area`, which lists them
north, west, south, east as the CDS expects.
"""
const AFRICA_BOX = (north = 40.0, west = -25.0, south = -36.0, east = 56.0)

"""
    forced(argv)

Return whether `--force` was given: rebuild or refetch an output that already exists.
"""
forced(argv) = "--force" in argv

"""
    repopath(path)

Return `path` relative to the repository root when it lies inside it, relative to the data
directory when it lies inside that, and otherwise its basename: a name that means the same on
every machine, for provenance records, which must carry no absolute path.
"""
function repopath(path::AbstractString)
    full = abspath(path)
    root = abspath(joinpath(@__DIR__, "..", ".."))
    startswith(full, root) && return relpath(full, root)
    data = datadir()
    startswith(full, data) && return joinpath("data", relpath(full, data))
    return basename(full)
end

# The script a provenance record is written by, named from the module `common.jl` was included
# into, so every sidecar says which file produced it without any script having to say so.
const SCRIPTS = Dict("ERA5Download" => "era5_download.jl",
                     "TwentyCRDownload" => "twentycr_download.jl",
                     "Inventory" => "inventory.jl",
                     "GBIFDownload" => "gbif_download.jl",
                     "GBIFClean" => "gbif_clean.jl",
                     "SpeciesTraits" => "species_traits.jl")
function _writer()
    return "data/src/" * get(SCRIPTS, string(nameof(@__MODULE__)),
               string(nameof(@__MODULE__)))
end

"""
    write_provenance(path; source, fields...)

Write `path.provenance.toml` beside `path`: the script that wrote it, the `source` it came
from, when, its size, and any further `fields` (strings, numbers, vectors or dictionaries of
those). Nothing here may name an absolute path or a person: pass paths through
[`repopath`](@ref) and keep credentials and addresses out of `fields`.
"""
function write_provenance(path::AbstractString; source::AbstractString,
                          fields...)
    record = Dict{String, Any}("writer" => _writer(),
                               "source" => source,
                               "written" => string(now(UTC)))
    isfile(path) && (record["bytes"] = filesize(path))
    for (k, v) in fields
        record[string(k)] = v
    end
    open(path * ".provenance.toml", "w") do io
        return TOML.print(io, record)
    end
    return nothing
end

"""
    provenance_path(path)

Return the path of the provenance sidecar `write_provenance` writes for `path`.
"""
provenance_path(path::AbstractString) = path * ".provenance.toml"
