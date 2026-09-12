# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Fetch the ERA5 monthly means the Africa runs use, one variable per decade file, as described in
# `era5_request.toml`, from the Climate Data Store through CDSAPI.jl. Resumable: a file already
# present is skipped, so a run interrupted by the CDS queue is simply started again.
#
#     julia --project data/src/era5_download.jl [--force] [--dry-run]
#
# Needs a `~/.cdsapirc` holding the CDS URL and key, as the CDS API instructions describe;
# CDSAPI.jl reads it (or `CDSAPI_URL` and `CDSAPI_KEY` from the environment). Writes under
# `$AFRICA_DATADIR/era5/` (default `data/era5/`), and a `.provenance.toml` beside each file with
# the request as sent and the CDS job that served it. To change what is fetched, edit
# `era5_request.toml`; a file whose request changed must be removed or refetched with `--force`,
# since presence is the only thing checked.

module ERA5Download

using CDSAPI
using Dates
using TOML

include("common.jl")

const SPEC = joinpath(@__DIR__, "era5_request.toml")
const POLL_SECONDS = 20.0

"""
    last_complete_year(today = today())

Return the last year every month of which is available at the CDS: the previous calendar year.
"""
last_complete_year(today::Date = Dates.today()) = Dates.year(today) - 1

"""
    decades(first, last)

Return the `(start, stop)` year pairs of the decade files from `first` to `last` inclusive,
the last pair cut at `last`.
"""
function decades(first::Integer, last::Integer)
    out = Tuple{Int, Int}[]
    start = first
    while start <= last
        stop = min(start + 9, last)
        push!(out, (start, stop))
        start = stop + 1
    end
    return out
end

"""
    build_request(spec, variable, years)

Return the request dictionary handed to the CDS for one `variable` (an entry of the spec's
`variables` table) over `years`: the spec's fixed `request` fields plus `variable` and `year`
as the CDS wants them, strings in vectors.
"""
function build_request(spec::AbstractDict, variable::AbstractDict, years)
    request = Dict{String, Any}(spec["request"])
    request["variable"] = [variable["cds"]]
    request["year"] = [string(y) for y in years]
    return request
end

"""
    jobs(spec, outdir, last; force = false)

Return the files still to fetch as `(variable, years, path)` named tuples: every variable in
`spec` for every decade from the spec's first year to `last`, skipping paths that exist unless
`force`.
"""
function jobs(spec::AbstractDict, outdir::AbstractString, last::Integer;
              force::Bool = false)
    out = NamedTuple{(:variable, :years, :path)}[]
    for variable in spec["variables"],
        (y0, y1) in decades(spec["decades"]["first"], last)
        path = joinpath(outdir, "era5_$(variable["short"])_$(y0)s.nc")
        (isfile(path) && !force) && continue
        push!(out, (variable = variable, years = y0:y1, path = path))
    end
    return out
end

"""
    retrieve(dataset, request, path; wait = POLL_SECONDS)

Fetch `request` from the CDS `dataset` into `path` through `CDSAPI.retrieve`, which submits
the job, polls it every `wait` seconds and downloads the result. Return the job's identifier,
or `nothing` if the client did not report one.
"""
function retrieve(dataset::AbstractString, request::AbstractDict,
                  path::AbstractString; wait::Real = POLL_SECONDS)
    job = CDSAPI.retrieve(dataset, request, path, wait = wait)
    return job isa AbstractDict ? get(job, "jobID", nothing) : nothing
end

# The spec, and the fetch. `retrieve(dataset, request, path)` defaults to the wrapper above and
# is taken as an argument so the loop can be exercised without the network.
function fetch(specfile::AbstractString, outdir::AbstractString;
               force::Bool = false, dryrun::Bool = false, retrieve = retrieve)
    spec = TOML.parsefile(specfile)
    dataset = spec["dataset"]
    mkpath(outdir)
    lastyear = last_complete_year()
    todo = jobs(spec, outdir, lastyear, force = force)
    println(length(todo), " file(s) to fetch into ", outdir, " (",
            length(spec["variables"]), " variables, decades from ",
            spec["decades"]["first"], " to ", lastyear, ").")
    if dryrun
        for job in todo
            println("  ", basename(job.path), ": ", job.variable["cds"], " ",
                    first(job.years), "-", last(job.years))
        end
        return todo
    end
    for job in todo
        request = build_request(spec, job.variable, job.years)
        println("Requesting ", basename(job.path), " ...")
        flush(stdout)
        # Fetch to a temporary name, so an interrupted transfer never looks complete.
        tmp = job.path * ".part"
        jobid = retrieve(dataset, request, tmp)
        mv(tmp, job.path, force = true)
        write_provenance(job.path,
                         source = "Copernicus Climate Data Store, " * dataset,
                         client = "CDSAPI.jl " * string(pkgversion(CDSAPI)),
                         job = string(something(jobid, "")),
                         request = request)
        println("  wrote ", basename(job.path), " (",
                round(filesize(job.path) / 2^20, digits = 1), " MiB)")
        flush(stdout)
    end
    return todo
end

function main(argv)
    fetch(SPEC, joinpath(datadir(), "era5"), force = forced(argv),
          dryrun = "--dry-run" in argv)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ERA5Download.main(ARGS)
end
