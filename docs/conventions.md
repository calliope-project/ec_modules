# Requirements and conventions

The tools in this project should focus on being easily re-used by others.
To ensure replicability without over-complication, here is a list of simple rules to follow.

## File structure

>[!warning]
>[`snakemake`](https://snakemake.readthedocs.io/en/stable/) has guidelines on how to structure workflows and wrappers on their websites.
>Unfortunately, their documentation is often very difficult to follow.
>We do our best to provide practical advice here.

### Modules

Please follow the [folder structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility) suggested by `snakemake` within your specific module folder. This will ensure others do not have issues when importing your workflow.

```txt
modules/
┗ example_module/
  ┣ config/
┃ ┃ ┣ example_config.yaml
┃ ┃ ┗ schema.yaml
  ┣ resources/ ---------| For static non-configurable files
  ┣ results/   ---------| ALWAYS place rule outputs here!
  ┣ workflow/
┃ ┃ ┣ envs/      ---------| conda env file(s)
┃ ┃ ┣ report/    ---------| optional report functionality (see below)
┃ ┃ ┣ rules/     ---------| snakemake rules
┃ ┃ ┣ scripts/   ---------| code
┃ ┃ ┗ Snakefile  ---------| please avoid other names (no `myworkflow.smk`)
  ┣ LICENSE
  ┗ README.md
```

Some suggestions:

- Make use of snakemake's [report](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html) function, which will allow users to easily follow your workflow, evaluate runtimes and even visualise plots. See [here](https://snakemake.github.io/resources/report.html) for an example of how this looks.
- Avoid placing large files in `/resources`. Instead, download them from `zenodo`.
- Use stable URLs (like DOIs) for downloads whenever possible!

### Wrappers

We follow the same contributing guidelines as the [`snakemake-wrappers`](https://snakemake-wrappers.readthedocs.io/en/stable/contributing.html) project.

```txt
wrappers/
┗ example_wrapper/
  ┣ tests/             ---------| Simple integration test
┃ ┃ ┣ Snakefile
  ┣ environment.yaml   ---------| conda env with static tool versions
  ┣ meta.yaml          ---------| obligatory wrapper description (see below)
  ┗ wrapper.py         ---------| passes inputs/params/output location to the wrapped tool
```

Some suggestions:

- `\meta.yaml` is an obligatory file used by snakemake to prepare inputs for your wrapper. See [here](https://snakemake-wrappers.readthedocs.io/en/stable/contributing.html#meta-yaml-file) for more info.
- Avoid developing code in the wrapper folder! These are meant to interface with external tools.

## File formats

Always try to use the following file formats.

- **Gridded data**: use [`xarray`](https://docs.xarray.dev/en/stable/) compatible formats. Read [here](https://help.marine.copernicus.eu/en/articles/8176692-how-to-choose-between-netcdf-and-zarr-format-using-the-toolbox) for comparisons.
  - [`zarr`](https://zarr.readthedocs.io/en/stable/): a modern format with excellent multi-threading and metadata support.
  - [`netCDF`](https://en.wikipedia.org/wiki/NetCDF): a widely supported format, but less efficient in certain situations.
- **Tabular data**: follow [tidy](https://vita.had.co.nz/papers/tidy-data.pdf) practices (columns are single variables, rows are observations).
  - `csv`
- **Configuration files**: a single file per module/wrapper.
  - [`yaml`](https://yaml.org/): human-friendly data specification.

## Metadata conventions

To future proof your workflow, try to follow these simple rules.

### Units

Unit conversion tools are not strictly necessary.
However, we recommend to use [`pint`](https://pint.readthedocs.io/en/stable/) if your data will be unit-sensitive.
This library automatically takes care of unit metadata and conversions, and is compatible with with [`xarray`](https://github.com/xarray-contrib/pint-xarray) and [`pandas`](https://github.com/hgrecco/pint-pandas).

### Country names

For existing countries, please use the [ISO 3166-1](https://en.wikipedia.org/wiki/ISO_3166-1) standard.
English short name, alpha-2, alpha-3 or numeric codes are all valid (e.g., "Grenada"::"GD"::"GRD"::308).
For historic regions, use [ISO 3166-3](https://en.wikipedia.org/wiki/ISO_3166-3) instead (e.g., "Czechoslovakia"::"CSHH"::200).

We recommend using the [`pycountry`](https://github.com/pycountry/pycountry) package, which follows [ISO 3166](https://en.wikipedia.org/wiki/ISO_3166).

>[!tip]
>To make data easier to parse, try to do the following:
>
>- Use alpha-2 or alpha-3 instead of English names.
>- Include the ISO 3166 numeric code. These are stable and non-repeating, even for former countries.

### Regions and subdivisions

Any official standard is valid: [NUTS](https://ec.europa.eu/eurostat/web/nuts), [GADM](https://gadm.org/) codes or [ISO 3166-2](https://en.wikipedia.org/wiki/ISO_3166-2) subdivisions.
The dataset version (`NUTS2024`, `GADM v1.3`, `ISO 3166-2:2013`) should be specified in either metadata or a separate column.

>[!warning]
> Specifying dataset version is extremely important! Subregions tend to be updated regularly, and additional processing might be needed in cases were two workflows use different versions.

### Timeseries

Please ensure that your timeseries are in [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601) UTC format **with no offset** in the form YYYY-MM-DD hh:mm:ss.

### Coordinate reference systems

Different coordinate systems are necessary [depending on the purpose of the spatial data](https://www.esri.com/arcgis-blog/products/arcgis-pro/mapping/gcs_vs_pcs/).
For simplicity, try to follow these standards:

- **Geographic data** (preserving *position* matters) should be in [EPSG 4326](https://epsg.io/4326). Most code should, in the end, return data in this format.
- **Projected data** (preserving *distance* or *area* matters) should use the representation that [best fits](https://learn.arcgis.com/en/projects/choose-the-right-projection/) the needs of the calculation. For European data, [EPSG:3035](https://epsg.io/3035) should fit most needs.
