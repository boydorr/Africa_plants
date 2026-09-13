# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Fetch the 20CRv3 monthly means described in `twentycr_request.toml` from NOAA PSL: one global
# netCDF file per variable, 1806 to 2015, over plain HTTPS with no credentials. Resumable: a file
# already present is skipped, so an interrupted run is simply started again.
#
#     julia --project data/src/twentycr_download.jl [--force] [--dry-run]
#
# Writes under `$AFRICA_DATADIR/twentycr/SI-MO/`, keeping PSL's directory and file names, and a
# `.provenance.toml` beside each file with the URL it came from. The files are global and are
# not cut here; the extension's climate step windows them.

module TwentyCRDownload

using Downloads
using TOML

include("common.jl")

const SPEC = joinpath(@__DIR__, "twentycr_request.toml")

"""
    jobs(spec, outdir; subsets = spec["subsets"], force = false)

Return the files still to fetch as `(subset, variable, url, path)` named tuples: every variable
in `spec` for each of `subsets`, skipping paths that exist unless `force`. The URL is the spec's
base, the variable's directory with the subset appended (PSL's layout: `2mSI-MO`, `accumsSI`),
and the file name.
"""
function jobs(spec::AbstractDict, outdir::AbstractString;
              subsets = spec["subsets"], force::Bool = false)
    out = NamedTuple{(:subset, :variable, :url, :path)}[]
    for subset in subsets, variable in spec["variables"]
        url = join([
                       spec["base"],
                       variable["directory"] * subset,
                       variable["file"]
                   ],
                   "/")
        path = joinpath(outdir, subset, variable["file"])
        (isfile(path) && !force) && continue
        push!(out,
              (subset = subset, variable = variable, url = url, path = path))
    end
    return out
end

# A downloader held to HTTP/1.1: on some networks a multi-hundred-megabyte HTTP/2 stream is cut
# off part way with `INTERNAL_ERROR`, and HTTP/1.1 transfers the same files without incident.
function _http11()
    downloader = Downloads.Downloader()
    downloader.easy_hook = (easy, _) -> Downloads.Curl.setopt(easy,
                                                              Downloads.Curl.CURLOPT_HTTP_VERSION,
                                                              Downloads.Curl.CURL_HTTP_VERSION_1_1)
    return downloader
end

"""
    download(url, path; attempts = 5)

Fetch `url` to `path` over HTTP/1.1, resuming from whatever `path` already holds with a byte
range on each retry, up to `attempts` transfers, the last failure rethrown. A server that
ignores the range and sends the whole file restarts the file rather than appending to it.
"""
function download(url::AbstractString, path::AbstractString;
                  attempts::Integer = 5)
    downloader = _http11()
    for attempt in 1:attempts
        have = isfile(path) ? filesize(path) : 0
        headers = have > 0 ? ["Range" => "bytes=$(have)-"] :
                  Pair{String, String}[]
        try
            open(path, "a") do io
                response = Downloads.request(url, output = io,
                                             headers = headers,
                                             downloader = downloader,
                                             throw = true)
                if have > 0 && response.status == 200
                    # The whole file came back: what was appended is the file again.
                    seekstart(io)
                    truncate(io, 0)
                    error("server ignored the byte range; restarting the file")
                end
                return response.status in (200, 206) ||
                       error("HTTP $(response.status) for $url")
            end
            return path
        catch e
            attempt == attempts && rethrow()
            println("  attempt ", attempt, " failed (",
                    first(sprint(showerror, e), 120), "); resuming from ",
                    isfile(path) ? filesize(path) : 0, " bytes")
            flush(stdout)
        end
    end
end

# The spec, and the fetch. `download(url, path)` defaults to the retrying HTTP/1.1 downloader
# above and is taken as an argument so the loop can be exercised without the network.
function fetch(specfile::AbstractString, outdir::AbstractString;
               subsets = nothing, force::Bool = false, dryrun::Bool = false,
               download = download)
    spec = TOML.parsefile(specfile)
    todo = jobs(spec, outdir, subsets = something(subsets, spec["subsets"]),
                force = force)
    println(length(todo), " file(s) to fetch into ", outdir, ".")
    if dryrun
        for job in todo
            println("  ", job.subset, "/", job.variable["file"], " <- ",
                    job.url)
        end
        return todo
    end
    for job in todo
        mkpath(dirname(job.path))
        println("Fetching ", job.subset, "/", job.variable["file"], " ...")
        flush(stdout)
        # Fetch to a temporary name, so an interrupted transfer never looks complete.
        tmp = job.path * ".part"
        download(job.url, tmp)
        mv(tmp, job.path, force = true)
        write_provenance(job.path,
                         source = "NOAA PSL, 20th Century Reanalysis V3 monthly means",
                         url = job.url, subset = job.subset,
                         variable = job.variable["short"],
                         role = job.variable["role"],
                         citation = "Slivinski et al. (2019) Q. J. R. Meteorol. Soc. 145:2876-2908, " *
                                    "doi:10.1002/qj.3598; data from NOAA PSL, https://psl.noaa.gov/data/gridded/data.20thC_ReanV3.html")
        println("  wrote ", basename(job.path), " (",
                round(filesize(job.path) / 2^20, digits = 1), " MiB)")
        flush(stdout)
    end
    return todo
end

function main(argv)
    spec = TOML.parsefile(SPEC)
    named = filter(a -> a in spec["subsets"], argv)
    fetch(SPEC, joinpath(datadir(), "twentycr"),
          subsets = isempty(named) ? nothing : named, force = forced(argv),
          dryrun = "--dry-run" in argv)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    TwentyCRDownload.main(ARGS)
end
