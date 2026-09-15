# NEWS

- v0.2.0
  - Changed
    - Rebuilt on EcoSISTEM 0.8 and Julia 1.11: the package now decides a projected Africa
      grid (`africa_area`) at the finest resolution the memory allows (`choose_cellsize`),
      following EcoSISTEM's HPC example, and is configured from the environment.
    - JuliaDB, JuliaDBMeta, IndexedTables, JLD and AxisArrays no longer work with Julia
      v1.11 and/or are no longer part of EcoSISTEM.
  - Added
    - Simple testing and workflow.
    - Data scripts under `data/src/` that fetch ERA5 monthly means from the Climate Data
      Store and a GBIF occurrence download for Africa, inventory the files, clean the records
      and derive the per-species climate table the run is parameterised from.
    - Cleaning drops records within 1 km of a country centroid (`data/src/centroids.jl`) and
      records whose locality mentions a greenhouse, and flags records on a quarter-degree
      lattice, within 1 km of a botanic garden (`data/src/gardens.jl`) or naming a garden in
      their locality, so that each use of the records chooses what to exclude.
    - `records_africa.csv` is the seeding dataset: every cleaned record inside Africa, flags kept.
    - `data/src/evi_download.jl` and `evi_mean.jl` fetch MODIS/Terra monthly EVI for 2001 to 2025
      and average it into one layer, the vegetation proxy for the collection-effort correction.
    - The ERA5 and 20CRv3 downloads go through EcoSISTEM (`era5requests`, `TwentyCR`,
      `fetchfiles`), which writes the provenance records.
- v0.1.0
  - The code behind case study 2 of the EcoSISTEM arXiv paper, on EcoSISTEM 0.1.
