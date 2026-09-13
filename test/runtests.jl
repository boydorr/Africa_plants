# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Smoke test: a small synthetic species pool on the real Africa grid at the coarsest
# resolution, run for a year. Reads WorldClim's annual mean temperature (cached after the first
# download) to position the grid, since a study area made only of synthetic layers cannot be
# placed by a named region.

module TestAfrica

using Test
using Africa
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources

ENV["AFRICA_SCALE"] = "small"

const TEMPERATURE = SourceSpec(WorldClim{BioClim}, :bio1)
const SUNLIGHT = UniformSpec(1000.0kJ / km^2 / month_mean_duration,
                             axis = SolarRadiation)
const NUMSPECIES = 20

@testset "Configuration" begin
    @test smallscale()
    @test isnothing(requested_cellsize())
    withenv("AFRICA_CELLSIZE" => "50km") do
        @test requested_cellsize() == 50.0km
    end
    withenv("AFRICA_CELLSIZE" => "7km") do
        @test_throws ErrorException requested_cellsize()
    end
end

# --- the data scripts, on fabricated inputs ---------------------------------------------------
#
# The scripts under `data/src/` fetch nothing here; their pure functions are exercised on files
# made up on the spot, so the join and the cleaning are checked before any real download exists.

const DATASRC = joinpath(@__DIR__, "..", "data", "src")
include(joinpath(DATASRC, "era5_download.jl"))
include(joinpath(DATASRC, "twentycr_download.jl"))
include(joinpath(DATASRC, "inventory.jl"))
include(joinpath(DATASRC, "gbif_download.jl"))
include(joinpath(DATASRC, "gbif_clean.jl"))
include(joinpath(DATASRC, "species_traits.jl"))

using NCDatasets
using CSV
using DataFrames
# Not plain `using Dates`: it would shadow Unitful's `year` and `day`.
using Dates: DateTime, Month
import Dates
using TOML
using Statistics

# Write a small ERA5-shaped file: `variable` over a 3 x 4 lattice for `months` consecutive months
# from January `year0`, valued `base + t` in month t at every cell, with the CDS's own axis names
# and a `units` attribute; `NaN` in the north-west corner throughout to stand for the sea.
function fake_era5(path, variable, units, year0, months; base = 0.0)
    NCDataset(path, "c") do ds
        defDim(ds, "valid_time", months)
        defDim(ds, "latitude", 3)
        defDim(ds, "longitude", 4)
        t = defVar(ds, "valid_time", Float64, ("valid_time",),
                   attrib = [
                       "units" => "days since 1900-01-01",
                       "calendar" => "proleptic_gregorian"
                   ])
        dates = [DateTime(year0, 1, 1) + Month(k - 1) for k in 1:months]
        t[:] = Dates.value.(dates .- DateTime(1900, 1, 1)) ./ 86_400_000
        defVar(ds, "latitude", Float64, ("latitude",))[:] = [10.0, 9.75, 9.5]
        defVar(ds, "longitude", Float64, ("longitude",))[:] = [
            20.0,
            20.25,
            20.5,
            20.75
        ]
        v = defVar(ds, variable, Float32,
                   ("longitude", "latitude", "valid_time"),
                   attrib = ["units" => units], fillvalue = Float32(-9999))
        for k in 1:months
            slab = fill(Float32(base + k), 4, 3)
            v[:, :, k] = slab
        end
        return v[1, 1, :] = fill(Float32(-9999), months)
    end
    return path
end

@testset "Data scripts: ERA5 request" begin
    @test ERA5Download.last_complete_year(Dates.Date(2026, 9, 11)) == 2025
    @test ERA5Download.decades(1940, 1965) ==
          [(1940, 1949), (1950, 1959), (1960, 1965)]
    @test ERA5Download.decades(1940, 1940) == [(1940, 1940)]
    spec = TOML.parsefile(ERA5Download.SPEC)
    @test spec["dataset"] == "reanalysis-era5-single-levels-monthly-means"
    @test [v["short"] for v in spec["variables"]] ==
          ["t2m", "tp", "swvl1", "ssr"]
    request = ERA5Download.build_request(spec, spec["variables"][2], 1950:1959)
    @test request["variable"] == ["total_precipitation"]
    @test request["year"] == string.(1950:1959)
    @test request["product_type"] == ["monthly_averaged_reanalysis"]
    @test request["area"] == [40, -25, -36, 56]
    @test request["data_format"] == "netcdf"
    @test length(request["month"]) == 12
    dir = mktempdir()
    todo = ERA5Download.jobs(spec, dir, 1955)
    @test length(todo) == 4 * 2
    @test basename(todo[1].path) == "era5_t2m_1940s.nc"
    @test todo[2].years == 1950:1955
    touch(todo[1].path)
    @test length(ERA5Download.jobs(spec, dir, 1955)) == 7
    @test length(ERA5Download.jobs(spec, dir, 1955, force = true)) == 8
    # The fetch loop without the network: a stand-in `retrieve` that writes the file.
    calls = String[]
    fake = (name, request, path) -> (
                                     push!(calls,
                                           name * " " * request["variable"][1]);
                                     write(path, "x"))
    expected = length(ERA5Download.jobs(spec, dir,
                                        ERA5Download.last_complete_year()))
    done = ERA5Download.fetch(ERA5Download.SPEC, dir, retrieve = fake)
    @test length(done) == expected == length(calls)
    @test all(startswith("reanalysis-era5-single-levels-monthly-means "), calls)
    @test all(isfile(j.path) for j in done)
    @test all(isfile(ERA5Download.provenance_path(j.path)) for j in done)
    # A second pass finds nothing to do, and a dry run lists without fetching.
    @test isempty(ERA5Download.fetch(ERA5Download.SPEC, dir, retrieve = fake))
    rm(done[1].path)
    @test length(ERA5Download.fetch(ERA5Download.SPEC, dir, dryrun = true)) == 1
    @test !isfile(done[1].path)
end

@testset "Data scripts: 20CRv3 request" begin
    spec = TOML.parsefile(TwentyCRDownload.SPEC)
    @test spec["subsets"] == ["SI-MO"]
    @test [v["short"] for v in spec["variables"]] ==
          ["air2m", "prate", "dswrf", "uswrf", "soilw"]
    dir = mktempdir()
    todo = TwentyCRDownload.jobs(spec, dir)
    @test length(todo) == 5
    @test todo[1].url ==
          "https://downloads.psl.noaa.gov/Datasets/20thC_ReanV3/Monthlies/2mSI-MO/air.2m.mon.mean.nc"
    @test todo[1].path == joinpath(dir, "SI-MO", "air.2m.mon.mean.nc")
    @test todo[end].url ==
          "https://downloads.psl.noaa.gov/Datasets/20thC_ReanV3/Monthlies/subsfcSI-MO/soilw.mon.mean.nc"
    @test length(TwentyCRDownload.jobs(spec, dir, subsets = ["SI-MO"])) == 5
    # The fetch loop without the network.
    urls = String[]
    fake = (url, path) -> (push!(urls, url); write(path, "x"))
    done = TwentyCRDownload.fetch(TwentyCRDownload.SPEC, dir,
                                  subsets = ["SI-MO"], download = fake)
    @test length(done) == 5 == length(urls)
    @test all(isfile(j.path) for j in done)
    @test all(isfile(TwentyCRDownload.provenance_path(j.path)) for j in done)
    rec = TOML.parsefile(TwentyCRDownload.provenance_path(done[1].path))
    @test rec["url"] == urls[1] && rec["subset"] == "SI-MO"
    @test isempty(TwentyCRDownload.fetch(TwentyCRDownload.SPEC, dir,
                                         subsets = ["SI-MO"], download = fake))
end

@testset "Data scripts: provenance" begin
    root = abspath(joinpath(@__DIR__, ".."))
    @test Inventory.repopath(joinpath(root, "data", "era5")) ==
          joinpath("data", "era5")
    @test Inventory.repopath(joinpath(root, "data", "src", "common.jl")) ==
          joinpath("data", "src", "common.jl")
    @test Inventory.repopath("/somewhere/else/file.nc") == "file.nc"
    dir = mktempdir()
    out = joinpath(dir, "thing.csv")
    write(out, "a\n")
    Inventory.write_provenance(out, source = "a test", n = 3)
    rec = TOML.parsefile(Inventory.provenance_path(out))
    @test rec["writer"] == "data/src/inventory.jl"
    @test rec["source"] == "a test" && rec["n"] == 3 && rec["bytes"] == 2
    @test !haskey(rec, "file")
    @test !occursin(homedir(), read(Inventory.provenance_path(out), String))
end

@testset "Data scripts: inventory" begin
    dir = mktempdir()
    fake_era5(joinpath(dir, "era5_t2m_1940s.nc"), "t2m", "K", 1940, 24,
              base = 270.0)
    table = Inventory.inventory([dir])
    @test nrow(table) == 1
    row = table[1, :]
    @test row.variable == "t2m"
    @test row.units == "K"
    @test row.time_name == "valid_time"
    @test row.time_count == 24
    @test row.time_first == DateTime(1940, 1, 1)
    @test row.time_last == DateTime(1941, 12, 1)
    @test row.time_even == false      # real month lengths differ
    @test row.lat_regular == true && row.lon_regular == true
    @test row.lat_step ≈ -0.25 && row.lon_step ≈ 0.25
    @test row.fillvalue == "-9999.0"
end

@testset "Data scripts: GBIF request and cleaning" begin
    withenv("GBIF_USER" => "someone") do
        json = GBIFDownload.predicate_json("someone@example.org", from = 1940)
        @test occursin("\"KINGDOM_KEY\", \"value\": \"6\"", json)
        @test occursin("POLYGON((-25.0 -36.0, 56.0 -36.0, 56.0 40.0, -25.0 40.0, -25.0 -36.0))",
                       json)
        @test occursin("SIMPLE_CSV", json)
        @test occursin("\"greaterThanOrEquals\", \"key\": \"YEAR\", \"value\": \"1940\"",
                       json)
        @test !occursin("lessThanOrEquals", json)
        early = GBIFDownload.predicate_json("someone@example.org", from = 1806,
                                            to = 1939)
        @test occursin("\"greaterThanOrEquals\", \"key\": \"YEAR\", \"value\": \"1806\"",
                       early)
        @test occursin("\"lessThanOrEquals\", \"key\": \"YEAR\", \"value\": \"1939\"",
                       early)
    end
    # The recorded filter carries the science and nothing about a person.
    filt = GBIFDownload.filter_json(from = 1806, to = 1939)
    @test occursin("\"lessThanOrEquals\", \"key\": \"YEAR\", \"value\": \"1939\"",
                   filt)
    @test !occursin("creator", filt) && !occursin("@", filt)
    @test occursin(filt,
                   GBIFDownload.predicate_json("x@y.org", from = 1806,
                                               to = 1939,
                                               creator = "someone"))
    @test GBIFDownload.REQUESTS.era5 == (from = 1940, to = nothing)
    @test GBIFDownload.REQUESTS.early == (from = 1806, to = 1939)
    s = GBIFDownload.parse_status("""{"key":"abc","doi":"10.15468/dl.xyz","status":"SUCCEEDED","downloadLink":"https://x/abc.zip"}""")
    @test s.status == "SUCCEEDED" && s.doi == "10.15468/dl.xyz" &&
          s.link == "https://x/abc.zip"
    @test GBIFDownload.parse_status("{}") ==
          (status = nothing, doi = nothing, link = nothing, license = nothing)
    @test GBIFDownload.parse_status("""{"license":"http://creativecommons.org/licenses/by-nc/4.0/legalcode"}""").license ==
          "http://creativecommons.org/licenses/by-nc/4.0/legalcode"

    @test GBIFClean.cellkey(0.0, 0.0) == (4500, 9000)
    @test GBIFClean.cellkey(0.01, 0.01) == GBIFClean.cellkey(0.0, 0.0)
    @test GBIFClean.cellkey(0.021, 0.0) != GBIFClean.cellkey(0.0, 0.0)

    header = join(["gbifID", "speciesKey", "species", "decimalLatitude",
                      "decimalLongitude",
                      "year", "month", "basisOfRecord",
                      "coordinateUncertaintyInMeters",
                      "issue"], '\t')
    rows = [
        "1\t100\tAcacia tortilis\t-1.5\t36.8\t1990\t3\tPRESERVED_SPECIMEN\t100\t",
        "2\t100\tAcacia tortilis\t-1.5\t36.8\t1991\t\tPRESERVED_SPECIMEN\t100\t",    # no month
        "3\t101\tFicus sur\t5.0\t10.0\t2001\t7\tLIVING_SPECIMEN\t\t",                # garden material
        "4\t101\tFicus sur\t5.0\t10.0\t2001\t7\tHUMAN_OBSERVATION\t50000\t",         # too uncertain
        "5\t\t\t5.0\t10.0\t2001\t7\tHUMAN_OBSERVATION\t\t",                          # no species
        "6\t102\tAloe vera\t-33.9\t18.4\t1850\t11\tHUMAN_OBSERVATION\t\t",
        "7\t102\tAloe vera\t-33.91\t18.41\t1851\t1\tOBSERVATION\t\t"]               # in an excluded cell
    dir = mktempdir()
    input = joinpath(dir, "in.csv")
    write(input, header * "\n" * join(rows, "\n") * "\n")
    output = joinpath(dir, "records.csv")
    excluded = Set([GBIFClean.cellkey(-33.91, 18.41)])
    result = GBIFClean.clean_records(input, output, excluded = excluded)
    @test result.seen == 7
    @test result.kept == 2
    @test result.bydecade == Dict(1990 => 1, 1850 => 1)
    kept = CSV.read(output, DataFrame)
    @test kept.gbifID == [1, 6]
    @test names(kept) == string.(GBIFClean.COLUMNS)
end

@testset "Data scripts: species traits join" begin
    dir = mktempdir()
    # Two decades of four variables; month t carries the value base + t so a window is checkable.
    for (v, base) in (("t2m", 270.0), ("tp", 0.0), ("ssr", 1000.0), ("swvl1", 0.0))
        fake_era5(joinpath(dir, "era5_$(v)_1940s.nc"), v, "x", 1940, 12,
                  base = base)
        fake_era5(joinpath(dir, "era5_$(v)_1941s.nc"), v, "x", 1941, 12,
                  base = base + 12)
    end
    stacks = Dict(v => SpeciesTraits.load_stack(dir, v)
                  for v in SpeciesTraits.VARIABLES)
    t2m = stacks["t2m"]
    @test length(t2m.months) == 24
    @test t2m.months[1] == (1940, 1) && t2m.months[end] == (1941, 12)
    @test size(t2m.data) == (4, 3, 24)

    # December 1941 at a cell with data: the window is months 13..24, values 283..294.
    @test SpeciesTraits.window_summary(t2m, 9.5, 20.75, 1941, 12,
                                       reducer = minimum) == 283.0
    @test SpeciesTraits.window_summary(t2m, 9.5, 20.75, 1941, 12) ==
          mean(283:294)
    # Before a full window exists, or off the series, there is no answer.
    @test isnan(SpeciesTraits.window_summary(t2m, 9.5, 20.75, 1940, 11))
    @test isnan(SpeciesTraits.window_summary(t2m, 9.5, 20.75, 1939, 12))
    @test isnan(SpeciesTraits.window_summary(t2m, 9.5, 20.75, 1942, 1))
    # The sea cell (north-west corner) has no data at all.
    @test isnan(SpeciesTraits.window_summary(t2m, 10.0, 20.0, 1941, 12))

    records = DataFrame(gbifID = 1:4, speciesKey = [1, 1, 2, 2],
                        species = ["A a", "A a", "B b", "B b"],
                        decimalLatitude = [9.5, 9.5, 9.75, 10.0],
                        decimalLongitude = [20.75, 20.75, 20.5, 20.0],
                        year = [1941, 1940, 1941, 1941], month = [12, 12, 6, 6])
    climate = SpeciesTraits.record_climate(records, stacks)
    @test climate.tmin ≈ [283.0, 271.0, 277.0, NaN] nans = true
    @test climate.tp ≈ [mean(13:24), mean(1:12), mean(7:18), NaN] nans = true
    table = SpeciesTraits.species_table(climate)
    @test nrow(table) == 1            # B has one record with climate, below the minimum
    @test table.species == ["A a"]
    @test table.numrecords == [2]
    @test table.tmin_mean ≈ [277.0]
    @test table.tmin_sd ≈ [std([283.0, 271.0])]
    @test table.swvl1_mean ≈ [mean([mean(13:24), mean(1:12)])]
end

@testset "Data scripts: the Africa outline" begin
    # Natural Earth's continent outline on WorldClim's own 10 arcminute lattice, first fetched
    # and then cached. Coastal Lagos sits in a cell that is mostly sea and is kept by the
    # one-cell tolerance; Madrid, Riyadh and the open Atlantic are out.
    places = DataFrame(name = ["Nairobi", "Antananarivo", "Cape Town", "Lagos",
                           "Sahara", "Madrid", "Riyadh", "mid-Atlantic"],
                       decimalLatitude = [-1.29, -18.9, -33.93, 6.45, 23.0,
                           40.4,
                           24.7, -10.0],
                       decimalLongitude = [36.82, 47.5, 18.42, 3.4, 10.0, -3.7,
                           46.7, -20.0])
    @test SpeciesTraits.in_africa(places) ==
          [true, true, true, true, true, false, false, false]
end

@testset "Choosing a grid" begin
    # No outline: the bounding box costs no download, and the smoke test wants speed.
    within = africa_within(outline = false)
    chosen = choose_cellsize(NUMSPECIES, regime = TEMPERATURE,
                             supply = SUNLIGHT, within = within)
    @test chosen.fits
    @test chosen.cellsize == maximum(CANDIDATE_CELLSIZES)
    @test chosen.perspecies ==
          perspecies_bytes(chosen.cellsize, regime = TEMPERATURE,
                           supply = SUNLIGHT, within = within)

    area = africa_area(chosen.cellsize, regime = TEMPERATURE,
                       supply = SUNLIGHT, within = within,
                       verbosity = :silent)
    environment = GridHabitat(regime = TEMPERATURE, supply = SUNLIGHT,
                              area = area)
    # One Int64 per cell of the whole grid: the priced footprint is the built one.
    @test chosen.perspecies == sizeof(Int64) * length(environment.active)
    @test count(environment.active) > 0

    species = build_species(NUMSPECIES, tolerance = (295.0K, 50.0K),
                            toleranceaxis = Temperature,
                            demand = 10.0kJ / month_mean_duration,
                            demandaxis = SolarRadiation, dispersal = 150.0km,
                            abundance = 10^5, seed = 1)
    eco = build_ecosystem(species, environment, seed = 1)
    simulate!(eco, 1year - 1month_mean_duration, 1month_mean_duration)
    @test sum(eco.abundances.matrix) > 0
    @test count(>(0), sum(eco.abundances.matrix, dims = 2)) > 0
end

end
