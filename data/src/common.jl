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
    write_provenance(path; source, fields...)

Write `path.provenance.toml` beside `path`: who wrote it, the `source` it came from, when, its
size, and any further `fields` (strings, numbers, vectors or dictionaries of those).
"""
function write_provenance(path::AbstractString; source::AbstractString,
                          fields...)
    record = Dict{String, Any}("writer" => "Africa_plants data/src",
                               "source" => source,
                               "written" => string(now(UTC)),
                               "file" => basename(path))
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
