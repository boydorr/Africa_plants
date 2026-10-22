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
    @test africa_within() isa ShapeMaskSpec
    @test africa_within().rule isa AnyOverlap
    @test africa_within().shape isa NaturalEarthSpec
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
include(joinpath(DATASRC, "centroids.jl"))
include(joinpath(DATASRC, "gardens.jl"))
include(joinpath(DATASRC, "evi_download.jl"))
include(joinpath(DATASRC, "evi_mean.jl"))
include(joinpath(DATASRC, "species_traits.jl"))

using NCDatasets
using CSV
using DataFrames
# Not plain `using Dates`: it would shadow Unitful's `year` and `day`.
using Dates: DateTime, Month
import Dates
using TOML
using Statistics
using SparseArrays: sparse
using DimensionalData: DimensionalData

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
    @test ERA5Download.cdsarea() == (40, -25, -36, 56)
    dir = mktempdir()
    requests = ERA5Download.requests(dir, lastyear = 2025)
    @test length(requests) == length(ERA5Download.CODES) == 4
    @test [basename(r.path) for r in requests[1]] ==
          [["era5_t2m_$(d)s.nc" for d in 1940:10:2010]; "era5_t2m_2020-2025.nc"]
    tp = requests[2][2]
    @test tp.dataset == "reanalysis-era5-single-levels-monthly-means"
    @test tp.request["variable"] == ["total_precipitation"]
    @test tp.request["year"] == string.(1950:1959)
    @test tp.request["area"] == [40, -25, -36, 56]
    @test tp.request["product_type"] == ["monthly_averaged_reanalysis"]
    @test length(tp.request["month"]) == 12
    # A dry run reports presence without the network, since a CDS request states no size.
    touch(requests[1][1].path)
    specs = [SourceSpec(ERA, code, files = files)
             for (code, files) in zip(ERA5Download.CODES, requests)]
    dry = EcoSISTEM.fetchfiles(specs, dryrun = true)
    @test length(dry) == 4 * 9
    @test count(e -> e.present, dry) == 1
end

@testset "Data scripts: 20CRv3 request" begin
    dir = mktempdir()
    assets = TwentyCRDownload.assets(dir)
    @test length(assets) == 5
    @test assets[1].url ==
          "https://downloads.psl.noaa.gov/Datasets/20thC_ReanV3/Monthlies/2mSI-MO/air.2m.mon.mean.nc"
    @test assets[1].path == joinpath(dir, "SI-MO", "air.2m.mon.mean.nc")
    @test assets[end].url ==
          "https://downloads.psl.noaa.gov/Datasets/20thC_ReanV3/Monthlies/subsfcSI-MO/soilw.mon.mean.nc"
    # A file already present is recorded on first use and fetched from nowhere.
    mkpath(dirname(assets[1].path))
    write(assets[1].path, "not a netCDF file, and never opened")
    @test EcoSISTEM.assetpath(assets[1]) == assets[1].path
    record = provenance(assets[1].path)
    @test record.url == assets[1].url
    @test record.dataset == "TwentyCR" && record.code == "air"
    @test isnothing(record.fetched)
end

@testset "Data scripts: EVI" begin
    # Months come from the day of year, leap years included.
    @test EVIDownload.granule_month("MOD13C2.A2020001.061.2020328142931.hdf") ==
          Dates.Date(2020, 1, 1)
    @test EVIDownload.granule_month("MOD13C2.A2024061.061.2024100000000.hdf") ==
          Dates.Date(2024, 3, 1)
    @test EVIDownload.granule_month("MOD13C2.A2025060.061.2025100000000.hdf") ==
          Dates.Date(2025, 3, 1)
    @test_throws ErrorException EVIDownload.granule_month("MYD13C2.A2020001.061.2020328142931.hdf")
    # A month listed twice keeps its newer production.
    names = ["MOD13C2.A2025152.061.2025300000000.hdf",
        "MOD13C2.A2025121.061.2025170000000.hdf",
        "MOD13C2.A2025152.061.2025200000000.hdf"]
    @test EVIDownload.latest_per_month(names) ==
          ["MOD13C2.A2025121.061.2025170000000.hdf",
        "MOD13C2.A2025152.061.2025300000000.hdf"]
    base = "https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/MOD13C2.061/"
    table = DataFrame("Granule UR" => ["a", "b", "c", "d"],
                      "Online Access URLs" =>
                          [base * "x/MOD13C2.A2001001.061.2020061230446.hdf",
                              base * "y/MOD13C2.A2025152.061.2025200000000.hdf",
                              base *
                              "z/MOD13C2.A2025152.061.2025300000000.hdf," *
                              "https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-public/BROWSE.jpg",
                              base * "w/MOD13C2.A2026001.061.2026040000000.hdf"
                          ],
                      "Size" => [89.7, 95.0, 96.0, 90.0])
    granules = EVIDownload.parse_catalogue(table)
    @test [g.name for g in granules] ==
          ["MOD13C2.A2001001.061.2020061230446.hdf",
        "MOD13C2.A2025152.061.2025300000000.hdf"]
    @test granules[2].month == Dates.Date(2025, 6, 1)
    @test granules[2].megabytes == 96.0
    @test endswith(granules[2].url, "z/MOD13C2.A2025152.061.2025300000000.hdf")
    withenv("EARTHDATA_TOKEN" => "  a-token\n") do
        @test EVIDownload.token() == "a-token"
    end
    if !Sys.iswindows()
        withenv("EARTHDATA_TOKEN" => nothing, "HOME" => mktempdir()) do
            @test_throws ErrorException EVIDownload.token()
        end
    end
    dir = mktempdir()
    good = joinpath(dir, "good.hdf")
    write(good, UInt8[0x0e, 0x03, 0x13, 0x01, 0x00])
    bad = joinpath(dir, "bad.hdf")
    write(bad, "<html>login</html>")
    @test EVIDownload._ishdf4(good) && !EVIDownload._ishdf4(bad)

    # Averaging keeps good (0) and marginal (1) months, and leaves out fill, snow or ice and cloud.
    total = zeros(Float64, 2, 2)
    count = zeros(Int16, 2, 2)
    EVIMean.accumulate!(total, count, Int16[5000 -3000; 2000 1000],
                        Int8[0 0; 2 1])
    EVIMean.accumulate!(total, count, Int16[3000 4000; 2500 -3000],
                        Int8[1 0; 3 0])
    @test count == Int16[2 1; 0 1]
    evi = EVIMean.finish(total, count)
    @test evi[1, 1] ≈ 0.4f0 && evi[1, 2] ≈ 0.4f0 && evi[2, 2] ≈ 0.1f0
    @test isnan(evi[2, 1])
    @test_throws ErrorException EVIMean.accumulate!(total, count,
                                                    zeros(Int16, 3, 2),
                                                    zeros(Int8, 3, 2))
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
    Inventory.write_sidecar(out, source = "a test", n = 3)
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

    # Whole, half and quarter degrees are lattice points; quarter-degree square centres are not.
    @test GBIFClean.on_lattice(12.5, 30.0)
    @test GBIFClean.on_lattice(-12.25, 30.75)
    @test !GBIFClean.on_lattice(-33.875, 18.875)
    @test !GBIFClean.on_lattice(12.3, 30.0)
    @test !GBIFClean.on_lattice(-1.5, 36.8)

    places = GBIFClean.PlaceIndex([0.0, 0.0499, 39.0], [0.0, 0.0499, 10.0])
    @test length(places) == 3
    @test GBIFClean.near(places, 0.0, 0.005)             # 0.56 km
    @test !GBIFClean.near(places, 0.0, 0.02)             # 2.2 km
    @test GBIFClean.near(places, 0.0501, 0.0501)         # across a bucket edge
    @test GBIFClean.near(places, 39.0, 10.011)           # 0.95 km, longitude shrunk
    @test !GBIFClean.near(places, 39.0, 10.013)          # 1.12 km
    @test !GBIFClean.near(GBIFClean.PlaceIndex(Float64[], Float64[]), 0.0, 0.0)
    @test sprint(show, places) == "PlaceIndex(3 places, 1.0 km)"

    catalogue = joinpath(mktempdir(), "centroids.tsv")
    write(catalogue,
          join([
                   "iso3\tiso2\tgadm1\tgbif_name\ttype\tarea_sqkm\tdecimal_longitude\tdecimal_latitude\tsource_locality_name\tsource_reference\tsource",
                   "TUN\tTN\t\tTunisia\tPCLI\t163610\t9.5\t34.1\tTunisia\thttps://x\tGeoNames",
                   "NAM\tNA\t\tNamibia\tPCLI\t825615\t17.2\t-22.1\tNamibia\thttps://x\tTGN",
                   "KEN\tKE\tKEN.1_1\tBaringo\tADM1\t11015\t35.9\t0.6\tBaringo\thttps://x\tGeoNames"
               ],
               "\n") * "\n")
    institutions = DataFrame(name = ["Kirstenbosch National Botanical Garden",
                                 "Karoo National Botanical Garden",
                                 "University of Pretoria, Botany Department",
                                 "Jardin Botanique de Kisantu", "Lost garden"],
                             type = ["Botanic_garden", "Herbarium", "Herbarium",
                                 missing, "Botanic_garden"],
                             decimalLatitude = [
                                 -33.9875,
                                 -33.6,
                                 -25.75,
                                 -5.13,
                                 missing
                             ],
                             decimalLongitude = [
                                 18.43,
                                 19.46,
                                 28.23,
                                 15.1,
                                 missing
                             ],
                             country = ["ZAF", "ZAF", "ZAF", "COD", missing],
                             source = ["BCGI", "indexherbariorum",
                                 "indexherbariorum",
                                 "indexherbariorum", "BCGI"])
    gardens = Gardens.botanic_gardens(institutions)
    # Typed as a garden, or named as one; a botany department is neither, and a garden with no
    # coordinates cannot exclude anything.
    @test gardens.name == ["Kirstenbosch National Botanical Garden",
        "Karoo National Botanical Garden", "Jardin Botanique de Kisantu"]
    @test gardens.type == ["Botanic_garden", "Herbarium", ""]
    @test names(gardens) ==
          ["name", "type", "country", "source", "latitude", "longitude"]
    cc = Centroids.country_centroids(catalogue)
    @test cc.iso2 == ["TN", "NA"]
    @test cc.latitude == [34.1, -22.1] && cc.longitude == [9.5, 17.2]

    @test GBIFClean.greenhouse("Kew, in the Temperate House glasshouse")
    @test !GBIFClean.greenhouse(missing)
    @test GBIFClean.garden_named("Hort. Kirstenbosch. Rock Garden.")
    @test GBIFClean.garden_named("PRETORIA; BOTANIESE TUIN; KWEKERY")
    @test !GBIFClean.garden_named("Libreville, Arboretum de Sibang.")
    @test !GBIFClean.garden_named("Forest relicts and cultivated areas, Mundemba town.")
    @test !GBIFClean.garden_named(missing)

    header = join(["gbifID", "speciesKey", "species", "decimalLatitude",
                      "decimalLongitude", "year", "month", "basisOfRecord",
                      "coordinateUncertaintyInMeters", "locality", "issue"], '\t')
    rows = [
        "1\t100\tAcacia tortilis\t-1.5\t36.8\t1990\t3\tPRESERVED_SPECIMEN\t100\tNairobi\t",
        "2\t100\tAcacia tortilis\t-1.5\t36.8\t1991\t\tPRESERVED_SPECIMEN\t100\t\t",   # no month
        "3\t101\tFicus sur\t5.0\t10.0\t2001\t7\tLIVING_SPECIMEN\t\t\t",               # garden material
        "4\t101\tFicus sur\t5.0\t10.0\t2001\t7\tHUMAN_OBSERVATION\t50000\t\t",        # too uncertain
        "5\t\t\t5.0\t10.0\t2001\t7\tHUMAN_OBSERVATION\t\t\t",                         # no species
        "6\t102\tAloe vera\t-33.9\t18.4\t1850\t11\tHUMAN_OBSERVATION\t\t\t",
        "7\t102\tAloe vera\t-33.91\t18.41\t1851\t1\tOBSERVATION\t\t\t",               # at a country centroid
        "8\t103\tProtea cynaroides\t-34.0\t18.5\t1990\t2\tPRESERVED_SPECIMEN\t\t\t",  # on the lattice
        "9\t103\tProtea cynaroides\t-34.125\t18.375\t1995\t2\tPRESERVED_SPECIMEN\t\tHort. Kirstenbosch.\t",
        "10\t104\tCycas revoluta\t-25.7\t28.2\t1996\t2\tHUMAN_OBSERVATION\t\tPretoria, greenhouse\t"
    ]
    dir = mktempdir()
    input = joinpath(dir, "in.csv")
    write(input, header * "\n" * join(rows, "\n") * "\n")
    output = joinpath(dir, "records.csv")
    result = GBIFClean.clean_records(input, output,
                                     centroids = GBIFClean.PlaceIndex([-33.91],
                                                                      [18.41]),
                                     gardens = GBIFClean.PlaceIndex([-34.126],
                                                                    [18.376]))
    @test result.seen == 10
    @test result.kept == 4
    @test result.centroid == 1 && result.greenhouse == 1
    @test result.lattice == 1 && result.garden == 1 && result.named == 1
    @test result.bydecade == Dict(1990 => 3, 1850 => 1)
    kept = CSV.read(output, DataFrame)
    @test kept.gbifID == [1, 6, 8, 9]
    @test kept.on_lattice == [false, false, true, false]
    @test kept.near_garden == [false, false, false, true]
    @test kept.garden_named == [false, false, false, true]
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

    # Kinshasa is on land; a point in the Gulf of Guinea is off the land, and one in the South
    # Atlantic off the grid altogether.
    places = grid_cells(environment, [-4.32, -2.0, -50.0], [15.31, 3.0, -20.0])
    @test places[1] > 0 && places[2] == 0 && places[3] == 0
    @test parent(environment.active)[places[1]]
    @test grid_cells(environment, [-2.0], [3.0], active = false)[1] > 0
    # Snapping moves a position in an inactive cell to the nearest active cell within reach.
    snapped = grid_cells(environment, [-2.0], [3.0], snap = 50)[1]
    @test snapped > 0 && parent(environment.active)[snapped]
    @test grid_cells(environment, [-4.32], [15.31], snap = 3) == places[1:1]

    mask = Bool[0 0 0 0; 0 0 0 1; 1 0 0 0]
    # From (1, 1) the active (3, 1) is 2 cells away and (2, 4) 3.2: the nearer wins; out of reach, none.
    @test Africa._nearestactive(mask, CartesianIndex(1, 1), 3) ==
          LinearIndices(mask)[3, 1]
    @test Africa._nearestactive(mask, CartesianIndex(1, 1), 1) == 0
    # From (2, 2) both are about as near: (3, 1) at 2 against (2, 4) at 4 in squared distance.
    @test Africa._nearestactive(mask, CartesianIndex(2, 2), 2) ==
          LinearIndices(mask)[3, 1]

    species = build_species(NUMSPECIES, tolerance = (295.0K, 50.0K),
                            toleranceaxis = Temperature,
                            demand = 10.0kJ / month_mean_duration,
                            demandaxis = SolarRadiation, dispersal = 150.0km,
                            abundance = 10^5, seed = 1)
    eco = build_ecosystem(species, environment, seed = 1)
    simulate!(eco, 1year - 1month_mean_duration, 1month_mean_duration)
    @test sum(eco.abundances.matrix) > 0
    @test count(>(0), sum(eco.abundances.matrix, dims = 2)) > 0
    # The totals kept for every recorded year are the landscape summed each way.
    totals = abundance_totals(eco)
    @test totals.species == vec(sum(eco.abundances.matrix, dims = 2))
    @test totals.cells == vec(sum(eco.abundances.matrix, dims = 1))
    @test sum(totals.species) == sum(totals.cells) == population(eco).total
end

@testset "The species pool" begin
    envelopes = DataFrame(species = ["A a", "B b", "C c"],
                          genus = ["A", "B", "C"], family = ["F", "F", "F"],
                          order = ["O", "O", "O"],
                          tmin = [280.0, 290.0, 300.0],
                          tp = [0.001, 0.002, 0.003],
                          ssr = [1.0e7, 1.5e7, 2.0e7],
                          swvl1 = [0.1, 0.2, 0.3])
    widths = DataFrame(species = ["A a", "B b", "C c", "D d"],
                       tmin_n = [2, 100, 1, 50],
                       tmin_sd = [1.0, 3.0, missing, 2.0],
                       tp_n = [2, 100, 40, 50],
                       tp_sd = [0.0, 0.002, 0.001, 0.001])
    built = species_pool(envelopes, widths)
    @test built.unmatched == ["D d"]
    # C has one record of temperature and so no width of its own there, and is kept on the prior.
    @test built.nowidth == ["C c"]
    @test built.pool.species == ["A a", "B b", "C c"]
    # No group holds five well-recorded species, so the prior is the median over the well recorded:
    # 3 K from B alone, and 0.0015 m/day from B and C. A's two records barely count against the
    # weight of twelve, a width of zero included; B's hundred mostly keep their own.
    @test built.pool.tmin_sd ≈ [sqrt((2 * 1.0^2 + 12 * 3.0^2) / 14), 3.0, 3.0]
    @test built.pool.tp ≈ [1.0, 2.0, 3.0]
    @test built.pool.tp_sd ≈
          [sqrt(12 * 1.5^2 / 14), sqrt((100 * 2.0^2 + 12 * 1.5^2) / 112),
        sqrt((40 * 1.0^2 + 12 * 1.5^2) / 52)]
    # With no weight on the prior the widths are the records' own.
    @test species_pool(envelopes, widths, strength = 0).pool.tmin_sd ≈
          [1.0, 3.0, 3.0]

    # A species takes the median of its well-recorded relatives, itself left out, from the first
    # level with enough of them: species 1 from its genus, species 3 - whose genus holds only
    # itself - from the family, and species 5, with no relatives, from everyone.
    taxonomy = [("g", "f"), ("g", "f"), ("h", "f"), ("g", "f"), ("i", "j")]
    priors = width_priors([1.0, 2.0, 4.0, 3.0, 9.0], [40, 40, 40, 40, 40],
                          taxonomy, members = 2)
    @test priors.prior ≈ [2.5, 2.0, 2.0, 1.5, 3.0]
    @test priors.depth == [1, 1, 2, 1, 3]
    # A poorly recorded species informs nobody, and a group it cannot name is passed over.
    priors = width_priors([1.0, 2.0, 4.0, missing], [40, 40, 5, 1],
                          [("g", "f"), ("g", "f"), ("g", "f"), (missing, "f")],
                          members = 2)
    @test priors.prior ≈ [1.5, 1.5, 1.5, 1.5]
    @test priors.depth == [3, 3, 1, 2]
    @test shrink_width(0.0, 2, 3.0, 12) ≈ sqrt(12 * 9 / 14)
    @test shrink_width(2.0, 5, 3.0, 0) ≈ 2.0

    species = africa_species(built.pool, density = 4.0 / m^2, seed = 3)
    @test species.names == ["A a", "B b", "C c"]
    again = africa_species(built.pool, density = 4.0 / m^2, seed = 3)
    @test species.movement.kernels[1].dist == again.movement.kernels[1].dist
end

@testset "Density widths" begin
    @test AnnualMinimumTemperature <: Temperature
    @test AnnualMeanPrecipitation <: Precipitation
    # Each width is written as a literal, so it is held to the pool it was worked out from wherever
    # that data is to hand: `sqrt(2 pi)` times the geometric mean of the pool's tolerance widths,
    # which puts a species of typical breadth at a suitability of one at its optimum.
    widths_file = joinpath(pkgdir(Africa), "data", "species_widths.csv")
    envelopes_file = get(ENV, "AFRICA_ENVELOPES",
                         expanduser("~/git/plant_trait/data/global_plant_climate_envelopes.csv"))
    if isfile(widths_file) && isfile(envelopes_file)
        pool = species_pool(CSV.read(envelopes_file, DataFrame),
                            CSV.read(widths_file, DataFrame)).pool
        typical(widths) = sqrt(2pi) * exp(sum(log, widths) / length(widths))
        @test ustrip(K, EcoSISTEM.densitywidth(AnnualMinimumTemperature)) ≈
              typical(pool.tmin_sd) rtol=0.05
        @test ustrip(mm / Unitful.d,
                     EcoSISTEM.densitywidth(AnnualMeanPrecipitation)) ≈
              typical(pool.tp_sd) rtol=0.05
    else
        @test_skip "the species widths and envelopes are not on this machine"
    end
end

@testset "Climate timing" begin
    # A two-year monthly stack whose value is its month's ordinal, on real month starts.
    dates = [DateTime(1940, 1, 1) + Month(k - 1) for k in 1:24]
    values = reshape(repeat(1.0:24.0, inner = 4), 2, 2, 24) .* K
    stack = EcoSISTEM.ClimateRaster(ERA,
                                    DimensionalData.DimArray(values,
                                                             (DimensionalData.Y(1:2),
                                                              DimensionalData.X(1:2),
                                                              DimensionalData.Ti(dates))),
                                    "t2m")
    @test Africa._retime(stack, ForwardClimate()) === stack
    reversed = Africa._retime(stack, ReversedClimate())
    @test reversed.array[1, 1, :] == (24.0:-1.0:1.0) .* K
    @test DimensionalData.lookup(reversed.array, DimensionalData.Ti) ==
          (0:23) .* 1.0month_mean_duration
    held = Africa._retime(stack, HeldClimateYear(1941))
    @test held.array[2, 2, :] == (13.0:24.0) .* K
    @test DimensionalData.lookup(held.array, DimensionalData.Ti) ==
          (1:12) .* 1.0month_mean_duration
    @test_throws ErrorException Africa._retime(stack, HeldClimateYear(1942))
    # The least of the twelve months to each, the first eleven taking the first full twelve.
    least = Africa._summarised(stack, TrailingMinimum())
    @test least.array[1, 1, :] == vcat(fill(1.0, 12), 2.0:13.0) .* K
    average = Africa._summarised(stack, TrailingMean())
    @test average.array[2, 1, :] ≈ vcat(fill(6.5, 12), 7.5:18.5) .* K
    @test Africa._summarised(stack, MonthlyValue()) === stack
    # A reversed run passes the same summaries backwards, and a held year is one value throughout.
    back = Africa._laidout(stack, TrailingMinimum(), ReversedClimate())
    @test parent(back.array)[1, 1, :] ==
          reverse(parent(least.array)[1, 1, :])
    year = Africa._laidout(stack, TrailingMean(), HeldClimateYear(1941))
    @test all(year.array[1, 1, :] .≈ 18.5K)
    @test DimensionalData.lookup(year.array, DimensionalData.Ti) ==
          (1:12) .* 1.0month_mean_duration
    # A cell with no data, as the sea is, stays without it and does not stop the rest.
    holed = copy(values)
    holed[1, 2, :] .= NaN * K
    sea = EcoSISTEM.ClimateRaster(ERA,
                                  DimensionalData.DimArray(holed,
                                                           DimensionalData.dims(stack.array)),
                                  "t2m")
    drowned = Africa._summarised(sea, TrailingMinimum())
    @test all(isnan, drowned.array[1, 2, :])
    @test drowned.array[1, 1, :] == least.array[1, 1, :]
    # A shorter run of months is summarised over what there is.
    @test Africa._summarised(stack, TrailingMean(3)).array[1, 1, 1:4] ≈
          [2.0, 2.0, 2.0, 3.0] .* K
    # Reversing twice returns the original slices.
    twice = Africa._retime(Africa._retime(stack, ReversedClimate()),
                           ReversedClimate())
    @test parent(twice.array) == parent(stack.array)
end

@testset "Filling to capacity" begin
    counts = [2 0 1; 6 0 0]
    # Cell 1: space for 12 at a mean demand of 1 m^2, light for 25.6 at a mean of 0.25 x 1 + 0.75 x
    # 0.5 kJ/day, so space binds and the shares of 1:3 fill 12. Cell 3 holds species 1 only and
    # its light, at 1 kJ/day each, binds at 4.
    supply = ([12.0, 5.0, 100.0] .* m^2, [16.0, 5.0, 4.0] .* kJ / Unitful.d)
    demand = ([1.0, 1.0] .* m^2, [1.0, 0.5] .* kJ / Unitful.d)
    filled = fill_to_capacity(sparse(counts), supply, demand)
    @test filled.species == [1, 2, 1]
    @test filled.cell == [1, 1, 3]
    @test filled.count == [3, 9, 4]
    # The same fact in other units fills the same.
    @test fill_to_capacity(sparse(counts),
                           (supply[1], uconvert.(J / Unitful.d, supply[2])),
                           demand) == filled
end

@testset "Seeding scale" begin
    area = fill(1.0e6m^2, 2, 2)
    evi = [0.5 0.25; NaN 0.0]
    fraction = [1.0 0.5; 1.0 1.0]
    # Two cells carry vegetation: 0.5 x 1e6 and 0.25 x 0.5 x 1e6, together 0.625e6 m^2.
    c = seeding_scale(evi, fraction, area, 2.5e6m^2)
    @test c ≈ 4.0 / m^2
    @test seeding_scale(evi, fraction, area, 2.5e6m^2, density = 2.0 / m^2) ≈
          8.0 / m^2
    # A cell that alone held every record would seed the whole of its share of the land.
    records = [2 1; 3 0]
    factors = seeding_factors(c, evi, fraction, area, records)
    @test factors[1, 1] ≈ 1.0e6
    @test factors[1, 2] ≈ 0.5e6
    # No EVI, or none worth the name, seeds only the records themselves; no records, nothing.
    @test factors[2, 1] == 1.0
    @test isnan(factors[2, 2])
    @test sum(filter(!isnan, factors) .* [2, 3, 1]) ≈ 2.5e6 + 3
    @test seeding_factors(c, [1.0e-9;;], [1.0;;], [1.0m^2;;], [5;;]) == [1.0;;]
    @test_throws ErrorException seeding_scale([NaN;;], [1.0;;], [1.0m^2;;],
                                              1.0m^2)

    # A record in the burn-in's first month is at elapsed zero; one after it is before the run.
    @test months_back([2025, 2026, 1940], [12, 1, 1]) == [0, -1, 1031]
    # Without a floor the factors are the scale's own shares, and a cell with no EVI seeds nothing.
    raw = seeding_factors(c, evi, fraction, area, records, floor = 0.0)
    @test raw[1, 1] ≈ 1.0e6 && raw[2, 1] == 0.0 && isnan(raw[2, 2])

    # One scale takes relative factors to a chosen total, each still seeding at least one.
    scaled = scale_to_total([2.0, 1.0, 0.0], 10.0)
    @test sum(scaled.counts) ≈ 10.0
    @test scaled.counts ≈ [6.0, 3.0, 1.0]
    @test scale_to_total([2.0, 1.0], 1.5).counts == [1.0, 1.0]
    @test scale_to_total([2.0, 1.0], 1.5).scale == 0.0

    # A record counts fully at the start and halves for each half-life before it.
    @test recency_weight(0.0year) == 1.0
    @test recency_weight(-1.0year) == 1.0
    @test recency_weight(40.0year) ≈ 0.25
    @test recency_weight(10.0year, halflife = 10.0year) ≈ 0.5

    grid = (lon = [0.025, 0.075], lat = [0.075, 0.025])
    @test fine_cells([0.07, 0.01, 1.0], [0.01, 0.07, 0.0], grid) == [1, 4, 0]
end

end
