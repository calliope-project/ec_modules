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

- `results/` can hold snakemake [report](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html) files, which will allow users to easily follow your workflow, evaluate runtimes and even visualise plots. See [here](https://snakemake.github.io/resources/report.html) for an example of how this looks.
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

- **Units**: although not strictly necessary, we recommend to use [`pint`](https://pint.readthedocs.io/en/stable/) if your data will be unit-sensitive. This library automatically takes care of unit metadata and conversions, and is compatible with with [`xarray`](https://github.com/xarray-contrib/pint-xarray) and [`pandas`](https://github.com/hgrecco/pint-pandas).
- **Country names**: existing countries should use the [ISO 3166-1 alpha-3](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3) standard (e.g., `MEX`), while historic countries should follow [ISO 3166-3](https://en.wikipedia.org/wiki/ISO_3166-3) (e.g., `SUHH`). We recommend using the [`pycountry`](https://github.com/pycountry/pycountry) package.
- **Regions/subdivisions**: [NUTS](https://ec.europa.eu/eurostat/web/nuts) or [GADM](https://gadm.org/) ids should be used, and the dataset version (`NUTS2024`, `GADM v1.3`) should be specified in either metadata or a separate column.

>[!warning]
> Specifying dataset version is extremely important! Subregions tend to be updated regularly, and additional processing might be needed in cases were two workflows use different versions.
