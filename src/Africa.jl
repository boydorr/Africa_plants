# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Simulating the plant species of Africa with EcoSISTEM: case study 2 of the EcoSISTEM paper,
# rebuilt on the current package. The module holds what every run shares - the configuration
# read from the environment, the continent as a study area, and the choice of a resolution the
# machine can afford. The experiment scripts under `examples/` drive it.

module Africa

using EcoSISTEM
using RasterDataSources: EarthEnv, LandCover
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Rasters: EPSG
using MPI
using CSV
using DataFrames
using Distributions: Uniform
using Random
using Dates: Dates
using SparseArrays: SparseArrays
using JLD2: JLD2
using Diversity: Diversity, Metacommunity, norm_sub_alpha
using DimensionalData: DimensionalData
using NCDatasets: NCDatasets

# EcoSISTEM's memory guidance for large runs, shared with its HPC scripts rather than copied.
# It is an example file of that package, not its API, so the path is pinned by the EcoSISTEM
# commit this package depends on.
include(joinpath(pkgdir(EcoSISTEM), "examples", "HPC", "memory.jl"))

include("config.jl")
include("area.jl")
include("landscape.jl")
include("species.jl")
include("seeding.jl")
include("placement.jl")
include("burnin.jl")
include("run.jl")

export smallscale, requested_cellsize, africa_within, AFRICA_CRS,
       CANDIDATE_CELLSIZES,
       perspecies_bytes, choose_cellsize, africa_area,
       AVAILABLE_CLASSES, era5_dir, available_ground, africa_layers,
       AbstractClimateTiming, ForwardClimate, ReversedClimate, HeldClimateYear,
       ENVELOPES_DOI, species_pool, envelopes_record,
       africa_species,
       MIN_EVI_MONTHS, seeding_scale, seeding_factors, months_back, fine_cells,
       read_seeding_grid, seeding_rows, recency_weight, scale_to_total,
       grid_cells, fill_to_capacity, positioning_layers, africa_ecosystem,
       gbif_record,
       burnin_table, run_burnin!, load_abundances, abundance_table,
       continue_from!,
       DIVERSITY_ORDERS, population, metacommunity_diversity,
       record_diversity, record_diversity!

end
