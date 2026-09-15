# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Fetch the MODIS/Terra monthly vegetation index files (MOD13C2 v061, global at 0.05 degrees) that
# `evi_mean.jl` averages into the vegetation proxy for the collection-effort correction: one HDF4
# file per month from `FIRST_MONTH` to `LAST_MONTH`, from NASA's LP DAAC.
#
#     julia --project data/src/evi_download.jl [--dry-run]
#
# The list of files comes from NASA's public catalogue (CMR), which needs no login; the files need a
# NASA Earthdata token, read from `EARTHDATA_TOKEN` or `~/.earthdata_token` and sent as a bearer
# token. A file already present is not fetched again, an interrupted one resumes from its part file,
# and each gets a provenance sidecar with its URL and checksum. Writes under `$AFRICA_DATADIR/evi/`.

module EVIDownload

using CSV
using DataFrames
using Dates
using Downloads
using SHA

include("common.jl")

const SHORT_NAME = "MOD13C2"
const PRODUCT_VERSION = "061"
const FIRST_MONTH = Date(2001, 1, 1)
const LAST_MONTH = Date(2025, 12, 1)
const DOI = "10.5067/MODIS/MOD13C2.061"
const CITATION = "Didan, K. (2021). MODIS/Terra Vegetation Indices Monthly L3 Global 0.05Deg CMG V061 [Dataset]. NASA Land Processes Distributed Active Archive Center. https://doi.org/10.5067/MODIS/MOD13C2.061"
const CMR = "https://cmr.earthdata.nasa.gov/search/granules.csv"

# The four bytes every HDF4 file starts with.
const HDF4_SIGNATURE = UInt8[0x0e, 0x03, 0x13, 0x01]

"""
    granule_month(name)

Return the first day of the month a MOD13C2 file or granule `name` covers, read from its
`A<year><day of year>` field.
"""
function granule_month(name::AbstractString)
    m = match(r"MOD13C2\.A(\d{4})(\d{3})\.", name)
    isnothing(m) && error("not a MOD13C2 granule name: $name")
    day = Date(parse(Int, m[1])) + Day(parse(Int, m[2]) - 1)
    return Date(Dates.year(day), Dates.month(day), 1)
end

# The production timestamp of a MOD13C2 granule name, which orders two productions of one month.
function _production(name::AbstractString)
    m = match(r"\.061\.(\d{13})", name)
    isnothing(m) && error("no production timestamp in $name")
    return m[1]
end

"""
    latest_per_month(names)

Return one of `names` per month, the newest production where a month appears more than once, in
month order.
"""
function latest_per_month(names::AbstractVector{<:AbstractString})
    best = Dict{Date, String}()
    for n in names
        month = granule_month(n)
        if !haskey(best, month) || _production(n) > _production(best[month])
            best[month] = String(n)
        end
    end
    return [best[k] for k in sort(collect(keys(best)))]
end

"""
    parse_catalogue(table; first = FIRST_MONTH, last = LAST_MONTH)

Return the files to fetch from a CMR granule table (`granules.csv`, whose `Online Access URLs` column
lists each granule's links and `Size` its megabytes) as `(name, month, url, megabytes)` named
tuples, one per month from `first` to `last`, the newest production of a month listed twice, in
month order.
"""
function parse_catalogue(table::DataFrame; first::Date = FIRST_MONTH,
                         last::Date = LAST_MONTH)
    byname = Dict{String, Tuple{String, Float64}}()
    for (links, size) in zip(table[!, "Online Access URLs"], table[!, "Size"])
        urls = filter(u -> startswith(u, "https://") && endswith(u, ".hdf"),
                      strip.(split(coalesce(links, ""), ",")))
        isempty(urls) && continue
        url = String(urls[1])
        byname[basename(url)] = (url, Float64(coalesce(size, 0.0)))
    end
    names = filter(n -> first <= granule_month(n) <= last,
                   latest_per_month(collect(keys(byname))))
    return [(name = n, month = granule_month(n), url = byname[n][1],
             megabytes = byname[n][2]) for n in names]
end

# The catalogue's table of MOD13C2 granules from `first` to the end of `last`'s month.
function _cataloguetable(first::Date, last::Date)
    url = CMR * "?short_name=$SHORT_NAME&version=$PRODUCT_VERSION" *
          "&temporal=$(first)T00:00:00Z,$(Dates.lastdayofmonth(last))T23:59:59Z&page_size=2000"
    body = IOBuffer()
    Downloads.request(url, output = body, throw = true)
    return CSV.read(take!(body), DataFrame)
end

"""
    token()

Return the NASA Earthdata token, from `EARTHDATA_TOKEN` or else from `~/.earthdata_token`, erroring
with where to put one when neither holds it.
"""
function token()
    t = strip(get(ENV, "EARTHDATA_TOKEN", ""))
    if isempty(t)
        file = joinpath(homedir(), ".earthdata_token")
        isfile(file) && (t = strip(read(file, String)))
    end
    isempty(t) &&
        error("no NASA Earthdata token: set EARTHDATA_TOKEN, or write the token to " *
              "~/.earthdata_token readable only by you")
    return String(t)
end

# A downloader held to HTTP/1.1: on some networks a long HTTP/2 stream is cut off part way.
function _http11()
    downloader = Downloads.Downloader()
    downloader.easy_hook = (easy, _) -> Downloads.Curl.setopt(easy,
                                                              Downloads.Curl.CURLOPT_HTTP_VERSION,
                                                              Downloads.Curl.CURL_HTTP_VERSION_1_1)
    return downloader
end

# Whether the file at `path` starts with the HDF4 signature, so an error page saved in its place
# is caught.
function _ishdf4(path::AbstractString)
    return open(io -> read(io, 4), path) == HDF4_SIGNATURE
end

"""
    fetch_granule(url, path, token; attempts = 5)

Fetch `url` to `path` with the Earthdata `token` as a bearer token over HTTP/1.1, writing to
`path.part` and resuming from what it holds on each retry, up to `attempts` transfers, the last
failure rethrown. The file is renamed into place only once it is complete and is an HDF4 file; a
server that ignores the byte range restarts the part file.
"""
function fetch_granule(url::AbstractString, path::AbstractString,
                       token::AbstractString; attempts::Integer = 5)
    part = path * ".part"
    downloader = _http11()
    for attempt in 1:attempts
        have = isfile(part) ? filesize(part) : 0
        headers = ["Authorization" => "Bearer " * token]
        have > 0 && push!(headers, "Range" => "bytes=$(have)-")
        try
            open(part, "a") do io
                response = Downloads.request(url, output = io,
                                             headers = headers,
                                             downloader = downloader,
                                             throw = true)
                if have > 0 && response.status == 200
                    truncate(io, 0)
                    error("the server sent the whole file in answer to a byte range")
                end
                return response.status in (200, 206) ||
                       error("HTTP $(response.status) for $(basename(url))")
            end
            if !_ishdf4(part)
                rm(part)
                error("$(basename(url)) arrived as something other than an HDF4 file")
            end
            mv(part, path, force = true)
            return path
        catch e
            attempt == attempts && rethrow()
            println("\n  attempt ", attempt, " failed (",
                    first(sprint(showerror, e), 160), "); retrying")
            flush(stdout)
            sleep(5 * attempt)
        end
    end
end

# The provenance sidecar of one fetched granule.
function _record(path::AbstractString, granule)
    return write_provenance(path, source = "NASA LP DAAC, MOD13C2 v061",
                            url = granule.url, month = string(granule.month),
                            doi = DOI, citation = CITATION,
                            sha256 = bytes2hex(open(sha256, path)))
end

function main(argv)
    dir = joinpath(datadir(), "evi")
    mkpath(dir)
    granules = parse_catalogue(_cataloguetable(FIRST_MONTH, LAST_MONTH))
    todo = filter(g -> !isfile(joinpath(dir, g.name)), granules)
    gigabytes = round(sum((g.megabytes for g in todo), init = 0.0) / 1024,
                      digits = 1)
    println(length(granules), " months in the catalogue from ", FIRST_MONTH,
            " to ", LAST_MONTH, "; ", length(todo), " to fetch (", gigabytes,
            " GB) into ", dir)
    "--dry-run" in argv && return nothing
    for g in granules
        path = joinpath(dir, g.name)
        (isfile(path) && !isfile(provenance_path(path))) && _record(path, g)
    end
    tok = token()
    for (k, g) in enumerate(todo)
        path = joinpath(dir, g.name)
        print("  [", k, "/", length(todo), "] ", g.name, " ...")
        flush(stdout)
        fetch_granule(g.url, path, tok)
        _record(path, g)
        println(" ", round(filesize(path) / 2^20, digits = 1), " MiB")
        flush(stdout)
    end
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    EVIDownload.main(ARGS)
end
