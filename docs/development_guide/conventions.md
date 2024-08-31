# Conventions

The tools in this project should focus on being easily re-used by others.
We want to increase compatibility without over-complicating things, so here is a list of simple pragmatic rules to follow.

## File formats

To ensure modules, wrappers and user workflows interact as seamlessly as possible, we recommend the following.

### Configuration data

- We recommend to use [`.yaml`](https://yaml.org/) in combination with `snakemake`'s [schema validation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#validation).

### Geospatial data

- For processing geospatial data, we recommend using [`rasterio`](https://github.com/rasterio/rasterio), [`geopandas`](https://geopandas.org/en/stable/), and [`gregor`](https://github.com/jnnr/gregor) libraries.
- For **polygon** data: use the `.geojson` format when possible.
- For **raster** data: use the geo-referenced format `.geotiff`.

### Gridded data

- We recommend using [`xarray`](https://docs.xarray.dev/en/stable/) to handle this type of data.
- Save in [`netCDF`](https://en.wikipedia.org/wiki/NetCDF) format. It is widely supported by most libraries.
- Add the following metadata to each variable. This will avoid ambiguity, and `xarray` will [automatically use them for plotting](https://docs.xarray.dev/en/stable/getting-started-guide/quick-overview.html#attributes).
    - `units` (e.g., `"W m-3"`, or `1` if unitless).
    - Optionally, `long_name` (e.g., `"Solar Irradiance"`)
- We encourage developers to follow the [CF Metadata Conventions](https://cfconventions.org/cf-conventions/cf-conventions.html), but this is not obligatory.

### Tabular data

- Please follow [tidy](https://vita.had.co.nz/papers/tidy-data.pdf) practices (columns are single variables, rows are observations).
- We recommend [`pandas`](https://pandas.pydata.org/docs/) for tabular data.
- Always save tabular data in `.csv` (comma separated) format.
- Always add a `units` column if applicable.

## Metadata conventions

To future proof our workflows we follow a few simple rules. In brief:

### All data

- :white_check_mark: Use `snake_case` for filenames, configuration keys, columns or data attributes.
- :white_check_mark: Add a `units` column or attribute.
- :ballot_box_with_check: Optionally, add a `long_name` column.

### National and subnational data

- :white_check_mark: Add a `country_id` column or attribute in [alpha-3](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3) format.
- :white_check_mark: Add a `subregion_id` column or attribute with the subregion code.
- :white_check_mark: Add a `subregion_spec` column or attribute. This specifies the version of the subregion standard (e.g., `NUTS2024`, `GADM v1.3`, `ISO 3166-2:2013`), which is necessary since subregion codes [change quite often](https://ec.europa.eu/eurostat/web/nuts/history).

### Timeseries data

- :white_check_mark: Ensure your timeseries follows [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601) UTC format (e.g., 2024-08-01T15:00:00+00:00).

### Geographic data

- :white_check_mark: Use `x` | `y` instead of `longitude` | `latitude` (or similars) if possible.
- :white_check_mark: **Geographic data** (preserving *position* matters) should be in [EPSG:4326](https://epsg.io/4326).
- :white_check_mark: **Projected data** (preserving *distance* or *area* matters) should use the representation that best fits the needs of the calculation. For European data, [EPSG:3035](https://epsg.io/3035) should fit most needs. Otherwise, we recommend using the [EPSG Dataset](https://epsg.org/home.html) to choose the best fitting reference system.

## Code conventions

We mostly follow python's [PEP 8 style](https://peps.python.org/pep-0008/).
To make things easier, we leverage several tools that help 'clean up' the code automatically.

- [`ruff`](https://docs.astral.sh/ruff/) is our main `python` linter and formatter. The tool should automatically follow our configuration in the `pyproject.toml` file. You can run it with the following:

    ```shell
    # Check the code for smells and antipatterns
    ruff check --fix
    # Automatically format the code to fit our style
    ruff format
    ```

    !!! warning "Docstrings"
        We have setup `ruff` to enforce ['Google style`](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) docstrings.
        This is to ensure our code remains understandable over time.
        Please make sure to document your code!
        This can be easily achieved with the [autoDocstring extension for VSCode](https://github.com/NilsJPWerner/autoDocstring), if you use that IDE.

- [`snakefmt`](https://github.com/snakemake/snakefmt) is our formatter for `snakemake` code and files. You can run it with the following:

    ```shell
    # Format all files in a folder
    snakefmt modules/my-module/
    ```

- We leverage [pre-commit ci](https://pre-commit.ci/) for continuous integration. This will execute before each PR and automatically run `ruff`, `snakefmt`, fix trailing spaces and other things. Nevertheless, we recommend you to check your code with `ruff` regularly!
