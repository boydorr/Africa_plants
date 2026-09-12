# Data acquisition and preparation

Scripts that fetch and prepare every input the Africa runs need. They write under `data/`
(gitignored, or the directory named by `AFRICA_DATADIR`), and beside each output a
`*.provenance.toml` saying where it came from, when, and how big it is. Each script is
idempotent: a file already present is not fetched or rebuilt unless `--force` is given.

Run them from the repository root in this order. Nothing here is run by the test suite except
the pure functions, which `test/runtests.jl` exercises on fabricated inputs.

| step | script | needs | writes |
| --- | --- | --- | --- |
| 1 | `era5_download.jl` | a CDS account: `~/.cdsapirc`, read by CDSAPI.jl | `era5/era5_<var>_<decade>.nc`, one variable per decade file |
| 2 | `inventory.jl` | the netCDF files | `inventory.csv`: every file's variables, unit strings, time axis, lattice |
| 3 | `gbif_download.jl` | a GBIF account: `GBIF_USER`, `GBIF_PWD`, `GBIF_EMAIL` in the environment, or `user`, `pwd`, `email` in `~/.gbifrc` | two downloads, `gbif/era5/<key>.csv` (1940 onward, the ERA5 years) and `gbif/early/<key>.csv` (1806-1939, the years only 20CRv3 reaches), each with its DOI in its provenance |
| 4 | `gbif_clean.jl` | the GBIF CSVs | `records.csv` from the `era5` download and `records_early.csv` from `early`: one row per usable record, cleaned |
| 5 | `species_traits.jl` | `records.csv`, the ERA5 files, a cached WorldClim layer to place the continent outline | `species.csv`, `records_africa.csv`, `record_climate.csv` |

```sh
julia --project data/src/era5_download.jl    # hours, queued at the CDS; resumable
julia --project data/src/inventory.jl
julia --project data/src/gbif_download.jl --request-only   # submit both; GBIF emails when ready
julia --project data/src/gbif_download.jl                  # later: wait if needed, then fetch
julia --project data/src/gbif_clean.jl       # every download present; or name one
julia --project data/src/species_traits.jl
```

## Credentials

Never in the repository. The CDS key lives in `~/.cdsapirc` as the `cdsapi` client expects; the
GBIF username, password and the email GBIF notifies come from the environment or from
`~/.gbifrc`, a TOML file with `user`, `pwd` and `email` lines. A script that
cannot find its credentials says which are missing and stops before touching the network.

## What each output is

- **`era5/`**: ERA5 monthly means on single levels, 0.25 degrees, over the Africa box, netCDF,
  one variable per file per decade: `t2m` (2 m temperature, K), `tp` (total precipitation,
  m/day as the monthly mean of the daily accumulation), `swvl1` (volumetric soil water, layer 1,
  m^3/m^3), `ssr` (surface net solar radiation, J/m^2/day). The request is
  `era5_request.toml`, which is also the fixture EcoSISTEM's download code is built against.
- **`inventory.csv`**: what the files actually say, recorded before anything trusts them:
  variable names, `units` attributes verbatim, time axis calendar, first and last dates,
  spacing and whether it is even, latitude and longitude ranges and spacing and whether the
  latitude lattice is regular, grid-mapping attributes, `scale_factor`, `add_offset`,
  `_FillValue`.
- **`gbif/era5/`, `gbif/early/`**: two GBIF occurrence downloads of Plantae with coordinates in
  the Africa box, 1940 onward and 1806 to 1939, as the simple tab-separated CSV GBIF produces,
  each with its request key and DOI. Two rather than one so that a paper built on the ERA5 years
  alone has exactly the dataset it used to cite.
- **`records.csv`** (from `era5`) and **`records_early.csv`** (from `early`): `gbifID,
  speciesKey, species, decimalLatitude, decimalLongitude, year, month` for every record with a species-level name, coordinates without a flagged geospatial
  issue, a month and a year, an accepted basis of record, coordinate uncertainty below the
  threshold, and not sharing a 0.02 degree cell with a botanic garden or country centroid when
  `data/gardens.csv` and `data/centroids.csv` are present.
- **`records_africa.csv`**: the records of `records.csv` that fall in an active cell of a
  geographic grid cut to Natural Earth's continent-level Africa (land inside the outline). The
  burn-in seeds from these; they are placed on the run's own grid at run time, since that grid
  is chosen then, and the effort raster is counted from them at the same point.
- **`record_climate.csv`**: `records_africa.csv` with the twelve-month climate summaries at each
  record's ERA5 cell and date: `tmin` (minimum monthly 2 m temperature, K), `tp`, `ssr`,
  `swvl1` (means), `NaN` where the window is not covered.
- **`species.csv`**: per species with at least two records carrying climate, the record count
  and the mean and standard deviation across records of those summaries, which the run uses as
  niche means and widths and as resource demands: `tmin_mean, tmin_sd, tp_mean, tp_sd,
  ssr_mean, swvl1_mean`.
