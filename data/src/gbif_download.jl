# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Request, wait for and fetch the GBIF occurrence downloads: every georeferenced plant record in
# the Africa box, as GBIF's simple tab-separated CSV, in **two** downloads so that each has a
# DOI of its own to cite - `era5`, 1940 onward, the years the ERA5 climate covers, and `early`,
# 1806 to 1939, the years only the NOAA 20CRv3 reanalysis reaches. A paper built on the ERA5 run
# alone cites the first. GBIF prepares each file on its side, which takes minutes to hours, so
# the request key is kept and a rerun resumes.
#
#     julia --project data/src/gbif_download.jl [era5|early] [--request-only] [--force] [--key KEY]
#
# `--request-only` submits the requests and stops; GBIF prepares the files and emails when they
# are ready, and a later run without the flag picks up the saved keys, waits if it must, and
# fetches. That is how to use a poor connection: request now, transfer when it is good.
#
# Needs `GBIF_USER`, `GBIF_PWD` and `GBIF_EMAIL` in the environment, or `user`, `pwd` and
# `email` in `~/.gbifrc` (a free gbif.org account; the email is where GBIF sends the completion
# notice). Writes `gbif/<name>/<key>.zip`, unpacks it to
# `gbif/<name>/<key>.csv`, and records the download's DOI - which the paper must cite - in the
# provenance sidecar. Only the Julia standard library is used for the HTTP, so this runs in the
# package environment with nothing extra.

module GBIFDownload

using Base64
using Downloads

include("common.jl")

const API = "https://api.gbif.org/v1/occurrence/download"
const PLANTAE = 6            # GBIF backbone taxon key for kingdom Plantae
const POLL_SECONDS = 60

"""
    REQUESTS

The two downloads, by name: the year range each covers, `to = nothing` meaning to the present.
`era5` matches the ERA5 climate (1940 onward); `early` is the century before it that only the
NOAA 20CRv3 reanalysis reaches (from its first year, 1806).
"""
const REQUESTS = (era5 = (from = 1940, to = nothing),
                  early = (from = 1806, to = 1939))

"""
    filter_json(; from, to = nothing, box = AFRICA_BOX)

Return the JSON of the download's filter alone: kingdom Plantae, coordinates present and
without a flagged geospatial issue, status present, year at or after `from` (and at or before
`to`, when given), within `box`. This is the scientific content of a request, and what the
provenance records; the geometry is the box as a WKT polygon, counter-clockwise as GBIF
requires.
"""
function filter_json(; from::Integer, to = nothing, box = AFRICA_BOX)
    w, e, s, n = box.west, box.east, box.south, box.north
    polygon = "POLYGON(($w $s, $e $s, $e $n, $w $n, $w $s))"
    upper = isnothing(to) ? "" :
            """{"type": "lessThanOrEquals", "key": "YEAR", "value": "$to"},\n      """
    return """
    {
      "type": "and",
      "predicates": [
        {"type": "equals", "key": "KINGDOM_KEY", "value": "$PLANTAE"},
        {"type": "equals", "key": "HAS_COORDINATE", "value": "true"},
        {"type": "equals", "key": "HAS_GEOSPATIAL_ISSUE", "value": "false"},
        {"type": "equals", "key": "OCCURRENCE_STATUS", "value": "PRESENT"},
        {"type": "greaterThanOrEquals", "key": "YEAR", "value": "$from"},
        $(upper){"type": "within", "geometry": "$polygon"}
      ]
    }"""
end

"""
    predicate_json(email; from, to = nothing, box = AFRICA_BOX, creator = ENV["GBIF_USER"])

Return the JSON body of a download request: [`filter_json`](@ref)'s filter, as the simple CSV
format, with `email` notified on completion and `creator` the GBIF username. Sent to GBIF and
never written to disk, since it names a person.
"""
function predicate_json(email::AbstractString; from::Integer, to = nothing,
                        box = AFRICA_BOX,
                        creator::AbstractString = get(ENV, "GBIF_USER", ""))
    return """
    {
      "creator": "$creator",
      "notificationAddresses": ["$email"],
      "sendNotification": true,
      "format": "SIMPLE_CSV",
      "predicate": $(filter_json(from = from, to = to, box = box))
    }
    """
end

# The credentials: `GBIF_USER`, `GBIF_PWD`, `GBIF_EMAIL` from the environment, or the same three
# as `user`, `pwd`, `email` in `~/.gbifrc` (TOML, like `~/.cdsapirc` for the CDS); an error names
# exactly which is missing.
function _credentials()
    rc = joinpath(homedir(), ".gbifrc")
    file = isfile(rc) ? TOML.parsefile(rc) : Dict{String, Any}()
    get_(env, key) = haskey(ENV, env) ? ENV[env] : get(file, key, nothing)
    cred = (user = get_("GBIF_USER", "user"), pwd = get_("GBIF_PWD", "pwd"),
            email = get_("GBIF_EMAIL", "email"))
    missing_ = [String(k) for (k, v) in pairs(cred) if isnothing(v)]
    isempty(missing_) ||
        error("GBIF credentials missing: " * join(missing_, ", ") *
              ". Set GBIF_USER, GBIF_PWD and GBIF_EMAIL in the environment, or write " *
              "user, pwd and email into ~/.gbifrc. A free account at gbif.org provides them.")
    return cred
end

function _authheader(cred)
    return "Authorization" =>
        "Basic " * base64encode(cred.user * ":" * cred.pwd)
end

# POST the request for the years `range` covers; GBIF answers with the download key as the
# plain response body.
function request_download(cred, range)
    body = IOBuffer()
    Downloads.request(API * "/request", method = "POST",
                      input = IOBuffer(predicate_json(cred.email,
                                                      from = range.from,
                                                      to = range.to,
                                                      creator = cred.user)),
                      output = body,
                      headers = ["Content-Type" => "application/json",
                          _authheader(cred)])
    key = strip(String(take!(body)))
    isempty(key) && error("GBIF returned no download key")
    return key
end

"""
    parse_status(json)

Return `(status, doi, link, license)` read out of GBIF's download-status JSON with pattern
matching, each `nothing` if absent; `license` is the most restrictive licence among the
datasets the download drew on, which is the licence the download carries. Status is one of `PREPARING`, `RUNNING`, `SUCCEEDED`, `KILLED`,
`CANCELLED`, `FAILED`, `FILE_ERASED`.
"""
function parse_status(json::AbstractString)
    field(k) = (m = match(Regex("\"$k\"\\s*:\\s*\"([^\"]*)\""), json);
                isnothing(m) ? nothing : String(m.captures[1]))
    return (status = field("status"), doi = field("doi"),
            link = field("downloadLink"), license = field("license"))
end

# GET the status of download `key`.
function status(key::AbstractString)
    body = IOBuffer()
    Downloads.request(API * "/" * key, output = body)
    return parse_status(String(take!(body)))
end

# Poll until the download has succeeded or failed.
function await(key::AbstractString)
    while true
        s = status(key)
        s.status in ("SUCCEEDED",) && return s
        s.status in ("KILLED", "CANCELLED", "FAILED", "FILE_ERASED") &&
            error("GBIF download $key ended as $(s.status)")
        println("  $key is $(s.status); waiting $(POLL_SECONDS)s")
        sleep(POLL_SECONDS)
    end
end

# Fetch the finished download to `gbif/<key>.zip`, unpack the CSV beside it, and write the
# provenance. The transfer goes to a temporary name so an interrupted one never looks complete.
function fetch(key::AbstractString, s, dir::AbstractString, cred, name, range)
    zip = joinpath(dir, key * ".zip")
    tmp = zip * ".part"
    Downloads.download(something(s.link, API * "/request/" * key), tmp)
    mv(tmp, zip, force = true)
    run(`unzip -o -q $zip -d $dir`)
    csv = joinpath(dir, key * ".csv")
    isfile(csv) || error("expected $(basename(csv)) inside the GBIF zip")
    provenance!(csv, key, s, cred, name, range)
    return csv
end

# The provenance sidecar of a fetched download: its DOI, which the paper must cite, the years
# it covers, and the request as sent.
function provenance!(csv::AbstractString, key::AbstractString, s, cred, name,
                     range)
    write_sidecar(csv, source = "GBIF occurrence download $key ($name)",
                  doi = something(s.doi, "unknown"),
                  license = something(s.license, "unknown"),
                  years = string(range.from, "-",
                                 isnothing(range.to) ? "present" : range.to),
                  citation = "GBIF.org (" * string(Dates.year(Dates.now())) *
                             ") GBIF Occurrence Download https://doi.org/" *
                             something(s.doi, "unknown"),
                  filter = filter_json(from = range.from, to = range.to))
    return nothing
end

# One named download: request it (or resume it, if `gbif/<name>/request.key` or `--key` names
# one), wait, fetch.
function download!(name::Symbol, cred; force::Bool = false, key = nothing,
                   requestonly::Bool = false)
    range = REQUESTS[name]
    dir = joinpath(datadir(), "gbif", String(name))
    mkpath(dir)
    keyfile = joinpath(dir, "request.key")
    if isnothing(key) && isfile(keyfile)
        key = String(strip(read(keyfile, String)))
    end
    if !isnothing(key) && isfile(joinpath(dir, key * ".csv")) && !force
        csv = joinpath(dir, key * ".csv")
        if !isfile(provenance_path(csv))
            provenance!(csv, key, status(key), cred, name, range)
            println(name, ": wrote the missing provenance for ", key, ".csv")
        end
        println(name, ": already have ", key,
                ".csv; nothing to do (use --force to refetch).")
        return nothing
    end
    if isnothing(key) || force
        key = request_download(cred, range)
        write(keyfile, key)
        println(name, ": requested GBIF download ", key, " (", range.from, "-",
                isnothing(range.to) ? "present" : range.to, ")")
    else
        println(name, ": resuming GBIF download ", key)
    end
    if requestonly
        st = status(key)
        println(name, ": ", key, " is ", something(st.status, "unknown"),
                "; run again without --request-only to fetch it")
        return nothing
    end
    s = await(key)
    csv = fetch(key, s, dir, cred, name, range)
    println(name, ": fetched ", csv, " (DOI ", something(s.doi, "unknown"), ")")
    return csv
end

# The downloads named on the command line (both by default), in order.
function main(argv)
    named = [Symbol(a)
             for a in argv
             if !startswith(a, "--") &&
        Symbol(a) in keys(REQUESTS)]
    names = isempty(named) ? collect(keys(REQUESTS)) : named
    i = findfirst(==("--key"), argv)
    key = !isnothing(i) && i < length(argv) ? argv[i + 1] : nothing
    (!isnothing(key) && length(names) != 1) &&
        error("--key applies to one download; name it: era5 or early")
    cred = _credentials()
    for name in names
        download!(name, cred, force = forced(argv), key = key,
                  requestonly = "--request-only" in argv)
    end
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    GBIFDownload.main(ARGS)
end
