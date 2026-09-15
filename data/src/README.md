# Data acquisition and preparation

Scripts that fetch and prepare every input the Africa runs need. They write under `data/`
(gitignored, or the directory named by `AFRICA_DATADIR`), and beside each output a
`*.provenance.toml` saying where it came from, when, and how big it is. EcoSISTEM writes those
records for the ERA5 and 20CRv3 downloads, which go through its `fetchfiles`; the other scripts
write their own. Each script is idempotent: a file already present is not fetched again, and a
table already built is rebuilt only with `--force`.

Run them from the repository root in this order. Nothing here is run by the test suite except
the pure functions, which `test/runtests.jl` exercises on fabricated inputs.

| step | script | needs | writes |
| --- | --- | --- | --- |
| 1 | `era5_download.jl` | a CDS account: `~/.cdsapirc`, read by CDSAPI.jl | `era5/era5_<var>_<decade>s.nc`, one variable per decade file; a decade not yet complete is `era5_<var>_<first>-<last>.nc` |
| 1b | `twentycr_download.jl` | nothing: NOAA PSL serves 20CRv3 over plain HTTPS | `twentycr/SI-MO/`, five global monthly files spanning 1806 to 2015 (2 m temperature, precipitation rate, downward and upward shortwave, soil moisture by layer), for the extension before 1940 |
| 1c | `evi_download.jl` | a NASA Earthdata token: `EARTHDATA_TOKEN` in the environment, or `~/.earthdata_token` | `evi/MOD13C2.A<year><day>.061.<production>.hdf`, the MODIS/Terra monthly vegetation index files for 2001 to 2025, one per month, global at 0.05 degrees |
| 1d | `evi_mean.jl` | the MOD13C2 files | `evi/evi_mean_2001_2025.nc`: the mean EVI of good and marginal months at each cell, the stand-in for how many plants a cell holds: collection effort is records over EVI, so a record in a greener or less sampled cell seeds more individuals |
| 2 | `inventory.jl` | the netCDF files | `inventory.csv`: every file's variables, unit strings, time axis, lattice |
| 3 | `gbif_download.jl` | a GBIF account: `GBIF_USER`, `GBIF_PWD`, `GBIF_EMAIL` in the environment, or `user`, `pwd`, `email` in `~/.gbifrc` | two downloads, `gbif/era5/<key>.csv` (1940 onward, the ERA5 years) and `gbif/early/<key>.csv` (1806-1939, the years only 20CRv3 reaches), each with its DOI in its provenance |
| 3b | `centroids.jl` | nothing: GBIF's Catalogue of Centroids on GitHub, at a pinned commit | `centroids.csv`: country centroids, which `gbif_clean.jl` drops records near |
| 3c | `gardens.jl` | nothing: CoordinateCleaner's institutions gazetteer on GitHub, at a pinned commit | `gardens.csv`: botanic gardens, which `gbif_clean.jl` flags records near |
| 4 | `gbif_clean.jl` | the GBIF CSVs | `records.csv` from the `era5` download and `records_early.csv` from `early`: one row per usable record, cleaned, with three flag columns |
| 5 | `species_traits.jl` | `records.csv`, a cached WorldClim layer to place the continent outline | `records_africa.csv`, the records the model seeds from |

```sh
julia --project data/src/era5_download.jl    # hours, queued at the CDS; resumable; --dry-run lists
julia --project data/src/twentycr_download.jl # about 2 GB, no credentials; extension only
julia --project data/src/evi_download.jl     # about 29 GB, an Earthdata token; resumable
julia --project data/src/evi_mean.jl
julia --project data/src/inventory.jl
julia --project data/src/gbif_download.jl --request-only   # submit both; GBIF emails when ready
julia --project data/src/gbif_download.jl                  # later: wait if needed, then fetch
julia --project data/src/centroids.jl
julia --project data/src/gardens.jl
julia --project data/src/gbif_clean.jl       # every download present; or name one
julia --project data/src/species_traits.jl
```

## Credentials

Never in the repository. The CDS key lives in `~/.cdsapirc` as the `cdsapi` client expects; the
GBIF username, password and the email GBIF notifies come from the environment or from
`~/.gbifrc`, a TOML file with `user`, `pwd` and `email` lines. The NASA Earthdata token for the EVI
files comes from `EARTHDATA_TOKEN` or from `~/.earthdata_token`, a file holding only the token and
readable only by its owner; tokens expire after about 60 days. A script that
cannot find its credentials says which are missing and stops before touching the network.

## What each output is

- **`era5/`**: ERA5 monthly means on single levels, 0.25 degrees, over the Africa box, netCDF,
  one variable per file per decade: `t2m` (2 m temperature, K), `tp` (total precipitation,
  m/day as the monthly mean of the daily accumulation), `swvl1` (volumetric soil water, layer 1,
  m^3/m^3), `ssr` (surface net solar radiation, J/m^2/day). The request each file came from,
  and the CDS job that served it, are in its provenance record.
- **`twentycr/SI-MO/`**: 20CRv3 monthly means from NOAA PSL, global, PSL's own file names,
  each file 1806-01 to 2015-12, fetched from the URLs in EcoSISTEM's `TwentyCR` catalogue. PSL's
  `SI` directory serves byte-identical files, so it is not fetched; the quality boundaries (1836
  and 1981) lie inside each file. Not cut to Africa here.
- **`evi/`**: MOD13C2 v061 (MODIS/Terra Vegetation Indices Monthly L3 Global 0.05Deg CMG,
  doi:10.5067/MODIS/MOD13C2.061), one HDF4 file per month from January 2001 to December 2025, with
  a provenance sidecar holding its URL and checksum; where the catalogue lists a month twice the
  newer production is kept. `evi_mean_2001_2025.nc` averages each cell's EVI over the months the
  product's pixel reliability rates good (0) or marginal (1), leaving out snow or ice, cloud,
  estimated values and fill, and records how many months went into each cell.
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
  speciesKey, species, decimalLatitude, decimalLongitude, year, month, on_lattice, near_garden,
  garden_named` for every record with a species-level name, coordinates, a month and a year, an
  accepted basis of record and coordinate uncertainty below the threshold, that is not within
  1 km of a country centroid in `data/centroids.csv` and whose locality does not mention a
  greenhouse, glasshouse, hothouse or shade house. The flags mark what one use may exclude and
  another keep: `on_lattice` (both coordinates whole multiples of a quarter degree, a position
  set to a coarse grid or rounded; quarter-degree square centres, the unit southern African
  herbaria georeference to, are not flagged), `near_garden` (within 1 km of a botanic garden in
  `data/gardens.csv`) and `garden_named` (the locality names a botanic garden, Hort., a nursery or
  cultivation). The provenance counts the records dropped and flagged.
- **`centroids.csv`**: `iso2, name, source, latitude, longitude` for the country-level entries
  of GBIF's Catalogue of Centroids at the commit named in `centroids.jl`. Regenerated rather
  than committed, since the catalogue states no licence.
- **`gardens.csv`**: `name, type, country, source, latitude, longitude` for every botanic garden
  in CoordinateCleaner's `institutions` gazetteer at the commit named in `gardens.jl`: the rows
  typed as gardens and those of other types (mostly herbaria) whose name says botanic garden or
  arboretum. Regenerated rather than committed, since CoordinateCleaner is GPL-3.
- **`records_africa.csv`**: the records of `records.csv` that fall in an active cell of a
  geographic grid cut to Natural Earth's continent-level Africa (land inside the outline), flags
  included. This is the seeding dataset: every flagged record is kept for seeding, since a garden
  plant grows where it is recorded and a lattice point is close enough to place a seed in a grid
  cell. They are placed on the run's own grid at run time, since that grid is chosen then, and
  the effort raster is counted from them at the same point.
