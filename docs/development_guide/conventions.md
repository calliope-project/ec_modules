# Conventions

The tools in this project should focus on being easily reused by others.
We want to increase compatibility without over-complicating things, so here is a list of simple pragmatic rules to follow.

## File formats

To ensure modules, wrappers and user workflows interact as seamlessly as possible, we require using the the following file formats.

### Configuration data

Use [`.yaml`](https://yaml.org/) files for `snakemake` configuration and schema validation.
Aim to make configuration files succinct and robust against user mistakes.

If a configuration value requires specific **units**, it must be specified in its name explicitly.

??? example "Clean naming in configuration file"

    ```yaml
    # Set this to 'False' to enable user inputs
    use_default_user_resources: True

    # Our recommended settings
    maximum_installable_mw_per_km2:
      pv_on_tilted_roofs: 160 # from [@Gagnon:2016][@Klauser:2016], i.e. 16% efficiency
      pv_on_flat_areas: 80 # from [@Gagnon:2016][@Klauser:2016][@Wirth:2017]
      onshore_wind: 8 # from [@EuropeanEnvironmentAgency:2009]
      offshore_wind: 15 # from [@EuropeanEnvironmentAgency:2009]
    roof_share: # from [@Trondle:2019]
      east: 0.1660
      north: 0.1817
      south: 0.1821
      west: 0.1681
      flat: 0.3020
    ```

### Geospatial data

There are two main types of spatially explicit data.

- For **polygon** data: use the `.geojson` format when possible.
- For **raster** data: use the geo-referenced format GeoTIFF (`.tiff`).

We recommend using [`rasterio`](https://github.com/rasterio/rasterio), [`geopandas`](https://geopandas.org/), and [`gregor`](https://github.com/jnnr/gregor) to handle these files.

### Gridded data

Use the `.nc` [netCDF4](https://en.wikipedia.org/wiki/NetCDF) format.
Each variable should specify **units** and (optionally) **long_name** attributes, following [CF Metadata Conventions](https://cfconventions.org/cf-conventions/cf-conventions.html).

We recommend using [`xarray`](https://docs.xarray.dev/en/stable/) to handle this type of data.

### Tabular data

Use `.csv` (comma separated) format.
Tabular data should follow [tidy data](https://vita.had.co.nz/papers/tidy-data.pdf) conventions if possible.
This can be omitted in cases were all data is in the same unit for different regions (like timeseries).

Additionally, all tabular datasets must include a second header indicating **units**, using `"No Unit"` for qualitative or unitless cases.
This avoids ambiguity for users, and is easier to clean than adding it between brackets or parenthesis.

???+ example "Tidy dataset of annual energy generation"

    | Year    | Region | Energy_Generation |
    |---------|--------|-------------------|
    | No Unit | No Unit| MW                |
    | 2020    | North  | 4500              |
    | 2020    | East   | 4800              |
    | 2021    | North  | 4600              |
    | 2021    | East   | 4900              |
    | 2022    | North  | 4700              |
    | 2022    | East   | 5000              |

??? example "Timeseries dataset"

    | Time                | AUT       | BEL       |
    |---------------------|-----------|-----------|
    | No Unit             | MW        | MW        |
    | 2016-01-01 00:00:00 | 0.1584614 | 0.1725854 |
    | 2016-01-01 01:00:00 | 0.1580427 | 0.1733172 |
    | 2016-01-01 02:00:00 | 0.1619893 | 0.1753312 |

We recommend processing tabular data with [`pandas`](https://pandas.pydata.org/docs/), validating it with [`pandera`](https://pandera.readthedocs.io/) and using [`pint-pandas`](https://pint-pandas.readthedocs.io) for unit conversion.

## Metadata conventions

To future proof our workflows we follow a few simple rules. In brief:

### All data

- :white_check_mark: We prefer `snake_case` for filenames, configuration keys, variable names, etc.
- :white_check_mark: Always specify `units` (e.g., `"W m-3"`, or `No Unit` if unitless).
See [File formats](conventions.md#file-formats) for specifics.

### National and subnational ID codes

- :white_check_mark: Add a `country_id` column or attribute in [alpha-3](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3) format.
- :white_check_mark: Add a `shape_id` column or attribute with unique subregion codes.
- :white_check_mark: Add a `shape_spec` column or attribute. This specifies the version of the subregion standard (e.g., `NUTS2024`, `GADM v1.3`, `ISO 3166-2:2013`), which is necessary since subregion codes [change quite often](https://ec.europa.eu/eurostat/web/nuts/history).

Make use of the [`pycountry`](https://github.com/pycountry/pycountry) library to easily translate country names.

??? example "Example: IDs in a `.geojson` file"

    | country_id   | shape_id     | shape_spec        | geometry     |
    |--------------|--------------|-------------------|--------------|
    | DEU          | DE11         | NUTS2024          | POLYGON(...) |
    | DEU          | DE14         | NUTS2024          | POLYGON(...) |
    | DEU          | DE40         | NUTS2024          | POLYGON(...) |
    | TZA          | TZ-01        | ISO 3166-2:2019   | POLYGON(...) |
    | TZA          | TZ-02        | ISO 3166-2:2019   | POLYGON(...) |
    | TZA          | TZ-19        | ISO 3166-2:2019   | POLYGON(...) |

### Timeseries data

- :white_check_mark: Ensure your timeseries follows [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601) UTC format (e.g., 2024-08-01T15:00:00+00:00).

### Geographic data

- :white_check_mark: Use `longitude` | `latitude` to express position (i.e., avoid ambiguous values like `x` | `y`).
- :white_check_mark: **Geographic data**, where preserving *position* matters, should be in [EPSG:4326](https://epsg.io/4326).
- :white_check_mark: **Projected data**, preserving *distance* or *area* matters, should use the representation that best fits the needs of the calculation. For European data, [EPSG:3035](https://epsg.io/3035) should fit most needs. Otherwise, we recommend using the [EPSG Dataset](https://epsg.org/home.html) to choose the best fitting reference system.

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
