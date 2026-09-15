# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The niche widths of the Africa species: for each species in `records_africa.csv`, the standard
# deviation across its records of the twelve-month climate at each record, one width per tolerance
# axis, each over the records that axis trusts.
#
#     julia --project -t auto data/src/species_widths.jl [--force]
#
# Temperature is the minimum monthly 2 m temperature over a record's month and the eleven before it
# (K), and precipitation the mean total precipitation over the same window (m/day, ERA5's own unit),
# the same summaries the published envelopes averaged into their means. Temperature leaves out
# records on the quarter-degree lattice; precipitation also leaves out records near or naming a
# botanic garden, whose plants may be watered. A width over fewer than two records is missing.
# Writes `species_widths.csv`.

module SpeciesWidths

using CSV
using DataFrames
using Statistics

include("species_traits.jl")
using .SpeciesTraits: load_stack, window_summary, datadir, forced,
                      write_sidecar

const MIN_RECORDS = 2

"""
    record_windows(records, t2m, tp)

Return `(tmin, tp)`, one value per row of `records`: the minimum of the `t2m` stack and the mean of
the `tp` stack over the twelve months ending in the record's month at its cell, `NaN` where the
window is not covered.
"""
function record_windows(records::DataFrame, t2m, tp)
    n = nrow(records)
    tmin = fill(NaN, n)
    rain = fill(NaN, n)
    lat, lon = records.decimalLatitude, records.decimalLongitude
    yr, mo = records.year, records.month
    Threads.@threads for r in 1:n
        tmin[r] = window_summary(t2m, lat[r], lon[r], yr[r], mo[r],
                                 reducer = minimum)
        rain[r] = window_summary(tp, lat[r], lon[r], yr[r], mo[r])
    end
    return (tmin = tmin, tp = rain)
end

# The standard deviation of the non-`NaN` values of `v` where `keep` holds, with how many went in;
# `missing` below `MIN_RECORDS`.
function _width(v::AbstractVector{Float64}, keep::AbstractVector{Bool})
    good = [x for (x, k) in zip(v, keep) if k && !isnan(x)]
    return (n = length(good),
            sd = length(good) < MIN_RECORDS ? missing : std(good))
end

"""
    width_table(records)

Group `records`, which carry `species`, the flags `on_lattice`, `near_garden` and `garden_named`,
and the window summaries `tmin` and `tp`, by species name, and return one row per name: `species,
tmin_n, tmin_sd, tp_n, tp_sd`. A name is the unit the published envelopes are keyed on, so records
GBIF files under several keys of one name are pooled. Temperature widths leave out records on the
lattice; precipitation widths also leave out records near or naming a garden.
"""
function width_table(records::DataFrame)
    grouped = groupby(records, :species)
    table = combine(grouped) do g
        temperature = _width(g.tmin, .!g.on_lattice)
        rain = _width(g.tp, .!(g.on_lattice .| g.near_garden .| g.garden_named))
        return (tmin_n = temperature.n, tmin_sd = temperature.sd,
                tp_n = rain.n, tp_sd = rain.sd)
    end
    return sort!(table, :species)
end

function main(argv)
    out = joinpath(datadir(), "species_widths.csv")
    if isfile(out) && !forced(argv)
        println("Already have ", out,
                "; nothing to do (use --force to rebuild).")
        return nothing
    end
    records = CSV.read(joinpath(datadir(), "records_africa.csv"), DataFrame)
    filter!(r -> !ismissing(r.species), records)
    println("Read ", nrow(records), " named records")
    flush(stdout)
    era5 = joinpath(datadir(), "era5")
    windows = record_windows(records, load_stack(era5, "t2m"),
                             load_stack(era5, "tp"))
    records.tmin = windows.tmin
    records.tp = windows.tp
    table = width_table(records)
    CSV.write(out, table)
    write_sidecar(out,
                  source = "records_africa.csv joined to ERA5 monthly means",
                  species = nrow(table),
                  records_with_temperature = count(!isnan, windows.tmin),
                  records_with_precipitation = count(!isnan, windows.tp),
                  species_with_temperature_width = count(!ismissing,
                                                         table.tmin_sd),
                  species_with_precipitation_width = count(!ismissing,
                                                           table.tp_sd),
                  window_months = 12,
                  temperature = "minimum monthly t2m over the window, K; records on the quarter-degree lattice left out",
                  precipitation = "mean monthly tp over the window, m/day; records on the lattice, near a botanic garden or naming one left out")
    println(nrow(table), " species -> ", out)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    SpeciesWidths.main(ARGS)
end
