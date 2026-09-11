# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Simulating the plant species of Africa with EcoSISTEM: case study 2 of the EcoSISTEM paper,
# rebuilt on the current package. The module holds what every run shares - the configuration
# read from the environment, the continent as a study area, and the choice of a resolution the
# machine can afford. The experiment scripts under `examples/` drive it.

module Africa

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Rasters: EPSG
using MPI

# EcoSISTEM's memory guidance for large runs, shared with its HPC scripts rather than copied.
# It is an example file of that package, not its API, so the path is pinned by the EcoSISTEM
# commit this package depends on.
include(joinpath(pkgdir(EcoSISTEM), "examples", "HPC", "memory.jl"))


export smallscale, requested_cellsize, africa_within, AFRICA_CRS,
       CANDIDATE_CELLSIZES,
       perspecies_bytes, choose_cellsize, africa_area

end
