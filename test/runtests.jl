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
