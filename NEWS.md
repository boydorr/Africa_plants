# NEWS

- v0.2.0
  - Changed
    - Rebuilt on EcoSISTEM 0.8 and Julia 1.11: the package now decides a projected Africa
      grid (`africa_area`) at the finest resolution the memory allows (`choose_cellsize`),
      following EcoSISTEM's HPC example, and is configured from the environment.
    - JuliaDB, JuliaDBMeta, IndexedTables, JLD and AxisArrays no longer work with Julia
      v1.11 and/or are no longer part of EcoSISTEM.
    - The data scripts' sidecar writer is `write_sidecar`, since EcoSISTEM exports
      `write_provenance`.
  - Removed
    - The JuliaDB-era example scripts (`prepData.jl`, `climateData.jl`, `burnin.jl`,
      `runAfrica.jl`, `plots.jl`), which no longer ran; the `africa_*` scripts replace them.
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
    - `species_pool` and `africa_species` build the simulated species from the published
      climate envelopes and the niche widths `data/src/species_widths.jl` measures, each species
      demanding ground and sunlight scaled by one `density`; water limits species through the
      precipitation tolerance rather than as a resource.
    - `seeding_scale`, `seeding_factors` and `seeding_rows` turn occurrence records into the
      individuals a backwards burn-in adds, scaled by EVI on the fine grid
      `data/src/seeding_grid.jl` writes; `fill_to_capacity` sets the starting state, with shares
      weighted by EVI and `recency_weight`, and `scale_to_total` makes the burn-in's additions
      equal it across the continent.
    - `grid_cells` places latitude and longitude positions on the projected grid, snapping those in
      inactive cells to the nearest active one; runs simulate every cell Africa's outline touches,
      so coastal records keep their cells.
    - `africa_layers` builds the ERA5 climate and the two supplies forward, reversed or with one
      year held; `africa_ecosystem`, `burnin_table`, `run_burnin!`, `continue_from!` and
      `record_diversity!` run the burn-in and the forward run, driven by
      `examples/africa_burnin.jl` and `examples/africa_forward.jl`, which
      `examples/africa_threads.jl` and `examples/africa_mpi.jl` run in turn on threads or across
      MPI ranks.
    - `examples/africa_plots.jl` draws species richness in the first and last years, the change
      between them and the departure from the neutral control, in its own `examples/Project.toml`
      environment so the package needs no plotting dependency.
- v0.1.0
  - The code behind case study 2 of the EcoSISTEM arXiv paper, on EcoSISTEM 0.1.
