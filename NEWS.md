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
- v0.1.0
  - The code behind case study 2 of the EcoSISTEM arXiv paper, on EcoSISTEM 0.1.
