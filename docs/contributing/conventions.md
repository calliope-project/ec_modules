# Requirements and conventions

The tools in this project should focus on being easily re-used by others.
To ensure compatibility without over-complicating things, here is a list of simple pragmatic rules to follow.

## File formats

To ensure modules, wrappers and user workflows interact as seamlessly as possible, we recommend the following.

### Gridded data

- We recommend using [`xarray`](https://docs.xarray.dev/en/stable/) to handle this type of data.
- Save in [`netCDF`](https://en.wikipedia.org/wiki/NetCDF) format. It is widely supported by most libraries.
- Add the following metadata to each variable. This will avoid ambiguity, and `xarray` will [automatically use them for plotting](https://docs.xarray.dev/en/stable/getting-started-guide/quick-overview.html#attributes).
    - `units` (e.g., `"W m-3"`, or `1` if unitless).
    - Optionally, `long_name` (e.g., `"Solar Irradiance"`)

We encourage developers to follow the [CF Metadata Conventions](https://cfconventions.org/cf-conventions/cf-conventions.html), but this is not obligatory.

### Tabular data

- Please follow [tidy](https://vita.had.co.nz/papers/tidy-data.pdf) practices (columns are single variables, rows are observations).
- We recommend [`pandas`](https://pandas.pydata.org/docs/) for tabular data.
- Always save tabular data in `.csv` (comma separated) format.
- Always add a `units` column if applicable.

### Configuration data

- We recommend to use [`.yaml`](https://yaml.org/) in combination with schema validation.

## Metadata conventions

To future proof our workflows we follow a few simple rules. In brief:

### All data

- :white_check_mark: Use `snake_case` for metadata like column names or attributes names.
- :white_check_mark: Add a `units` column or attribute.
- :ballot_box_with_check: Optionally, add a `long_name` column.

### National and subnational data

- :white_check_mark: Add a `country_id` column or attribute in [alpha-3](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3) format.
- :white_check_mark: Add a `subregion_id` column with the subregion code.
- :white_check_mark: Specify the version of the subregion standard (e.g., `NUTS2024`, `GADM v1.3`, `ISO 3166-2:2013`) in a separate column. This is necessary since subregion codes [change often](https://ec.europa.eu/eurostat/web/nuts/history).

### Timeseries data

- :white_check_mark: Ensure your timeseries follows [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601) UTC format (e.g., 2024-08-01T15:00:00+00:00).

### Geographic data

- :white_check_mark: Use `x` | `y` instead of `longitude` | `latitude` (or similars).
- :white_check_mark: **Geographic data** (preserving *position* matters) should be in [EPSG:4326](https://epsg.io/4326).
- :white_check_mark: **Projected data** (preserving *distance* or *area* matters) should use the representation that [best fits](https://learn.arcgis.com/en/projects/choose-the-right-projection/) the needs of the calculation. For European data, [EPSG:3035](https://epsg.io/3035) should fit most needs.
