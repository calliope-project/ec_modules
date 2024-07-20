# The easy energy system workflows and wrappers repository

This aims to be a collection of helper workflows and wrappers that aid in quickly developing Energy System (ESM) models, inspired by Snakemake's [workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/) and [wrapper repository](https://snakemake-wrappers.readthedocs.io/en/stable/).
Its main purpose is enabling researchers to share data workflows between studies to avoid the "bloat" problem in large energy system workflows.

## Setup

If you wish to test the whole repository, please run.

```shell
mamba env create -f environment.yaml
mamba activate ec_modules
```

>[!note]
> The `ec_modules` environment is mostly used for unit testing! See below for development cases.

If you wish to develop / test an individual module or wrapper, please install its specific environment:

```shell
mamba env create -f wrappers/name_here/environment.yaml
mamba activate name_here
```

## Requirements and conventions

Successful workflows are those that can be easily re-used by others.
To ensure this without over-complicating things, here is a list of simple rules to follow.

### File structure

Follow these simple rules to ensure compatibility. See [`snakemake`](https://snakemake.readthedocs.io/en/stable/) for more info on how to develop maintainable workflows.

#### Modules

Please follow the [folder structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility) suggested by `snakemake` within your specific module folder. This will ensure others do not have issues when importing your workflow.

```txt
modules/
┗ example_module/
  ┣ config/
┃ ┃ ┣ example_config.yaml
┃ ┃ ┗ schema.yaml
  ┣ resources/ ---------| Static files needed by your workflow
  ┣ results/   ---------| ALWAYS place rule outputs here!
  ┣ workflow/  ---------| Executed code lives here
┃ ┃ ┣ envs/      ---------| conda env files
┃ ┃ ┣ report/    ---------| figures
┃ ┃ ┣ rules/     ---------| snakemake rules
┃ ┃ ┣ scripts/   ---------| code
┃ ┃ ┗ Snakefile
  ┣ LICENSE
  ┗ README.md
```

#### Wrappers

We follow the same contributing guidelines as the [`snakemake-wrappers`](https://snakemake-wrappers.readthedocs.io/en/stable/contributing.html) project.

```txt
wrappers/
┗ example_wrapper/
  ┣ tests/
┃ ┃ ┣ Snakefile
  ┣ environment.yaml
  ┣ meta.yaml
  ┗ wrapper.py
```

### File format conventions

Always try to use the following file formats.

- **Gridded data**: use [`xarray`](https://docs.xarray.dev/en/stable/) compatible formats. Read [here](https://help.marine.copernicus.eu/en/articles/8176692-how-to-choose-between-netcdf-and-zarr-format-using-the-toolbox) for comparisons.
  - [`zarr`](https://zarr.readthedocs.io/en/stable/): a modern format with excellent multi-threading and metadata support.
  - [`netCDF`](https://en.wikipedia.org/wiki/NetCDF): a widely supported format, but less efficient in certain situations.
- **Tabular data**: follow [tidy](https://vita.had.co.nz/papers/tidy-data.pdf) practices (columns are single variables, rows are observations).
  - `csv`
- **Configuration files**: a single file per module/wrapper.
  - [`yaml`](https://yaml.org/): human-friendly data specification.

### Metadata conventions

To future proof your workflow, try to follow these simple rules.

- **Units**: although not strictly necessary, we recommend to use [`pint`](https://pint.readthedocs.io/en/stable/) if your data will be unit-sensitive. This library automatically takes care of unit metadata and conversions, and is compatible with with [`xarray`](https://github.com/xarray-contrib/pint-xarray) and [`pandas`](https://github.com/hgrecco/pint-pandas).
- **Country names**: existing countries should use the [ISO 3166-1 alpha-3](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3) standard (e.g., `MEX`), while historic countries should follow [ISO 3166-3](https://en.wikipedia.org/wiki/ISO_3166-3) (e.g., `SUHH`). We recommend using the [`pycountry`](https://github.com/pycountry/pycountry) package.
- **Regions/subdivisions**: [NUTS](https://ec.europa.eu/eurostat/web/nuts) or [GADM](https://gadm.org/) ids should be used, and the dataset version should be specified in either metadata or a separate column (e.g., `NUTS2024 | DE24B`).

>[!warning]
> Specifying dataset version is extremely important! Subregions tend to be updated regularly, and additional processing might be needed in cases were two workflows use different versions.

- **Outputs should only be files**: meaning that code variables cannot be passed between modules.
- **Outputs must be specified by templates**: modules might need to use the output of other modules as input (e.g., a _pv_capacity_ module might use a `*.geojson` file produced by a _regions_ module). To avoid breaking functionality, output template files must be specified in the `templates/` folder of a module to let other developers know what they'll be working with.
- **File structure should be standardised**: to help users reach relevant files quickly, modules must at least have _src_, _out_, _templates_ and _tmp_ folders.
- **Containerised**: each module must use its own _env_ to keep version conflicts at a minimum.
- **Isolated configuration YAML**: some sort of module configurability is needed to keep things flexible. However, modules must only have access to their own configuration parameters to avoid potential conflicts. Providing template configuration files is also recommended.
- **Debuggable and testable**: problems are bound to happen, so modules should be relatively easy to debug, and should include a few test cases by default.

![Wrapper example](docs/images/module_setup.png)


# Explanation

## What are energy model workflows?

In recent years, energy modelling researchers have made efforts to improve the transparency of how their models are created.
Multiple open-source model-building workflows have appeared, each delineating the many steps and tools employed to build a model.
These workflows usually aim to be compatible with a single modelling framework (e.g., a "language"), and are composed of a series of "steps", usually taking the form of a direct acyclic graph (DAG).

- [Sector-Coupled Euro Calliope](https://github.com/calliope-project/sector-coupled-euro-calliope): a workflow that builds a high-resolution model of Europe for the Calliope framework. Built in Snakemake, currently has 82 steps.
- [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur/tree/master): another european model workflow, built for the PyPSA framework. Also built in Snakemake, currently with 105 steps.
- [message-ix-models](https://docs.messageix.org/projects/models/en/latest/): a data processing toolbox to build MESSAGEix-GLOBIOM models. Number of steps unknown.
- [OSeMOSYS Global](https://github.com/OSeMOSYS/osemosys_global): workflow to build global electricity system models, built for the OSyMOSYS framework. Currently has 26 steps.
- [SpineToolbox](https://github.com/spine-tools/SpineOpt.jl): a toolbox to build SpineOpt models. GUI based.

![Example workflow](https://raw.githubusercontent.com/calliope-project/sector-coupled-euro-calliope/main/rulegraph.png)

## What is bloat?

In theory, more research should lead to scientific refinement and improved insights.
Research _needs_ previous work to be accessible in order to advance the field and produce better insights.
The current approach of model-specific workflows hampers progress in two ways.

1. By making sharing methodological improvements more difficult, since these processes are "locked-in" to a specific tool (e.g., Calliope, PyPSA, OSeMOSYS, etc). This leads to a lot of re-implementation between communities, which is error prone and more difficult to evaluate, which diminishes trust.
2. By making model workflows less comprehensible over time. Workflows they tend to grow in size and complexity as more studies are conducted with them, eventually turning into [black boxes](https://doi.org/10.1088/2516-1083/ad371e), increasing the risk of combining incompatible assumptions or using depreciated data, which in turn can lead to misleading studies.

![Bloat example](docs/images/bloat-growth-problem.png)

## Looking at other fields for a solution

Energy research is hardly the only field were this is a problem: bioinformatics and atmospheric modelling face similar problems.
Both of these fields share common [data characteristics](https://www.bsiranosian.com/bioinformatics/why-are-bioinformatics-workflows-different/) with energy research:

- Large file sizes.
- The need to process both common and proprietary formats.
- Relatively high computational intensity when processing data.
- Heterogeneous and complicated data processing steps.

The widespread use of Snakemake in energy model workflows hints at current issues being caused by poor collaboration, not by any inherent differences between these fields.
Snakemake was specifically developed to make data processing methods easier between bioinformaticians.
**Better workflow practices is needed to abstract complexity, letting people focus on doing science instead of reinventing the wheel!**

## Further reading:

- Wratten, L., Wilm, A. & Göke, J. Reproducible, scalable, and shareable analysis pipelines with bioinformatics workflow managers. Nat Methods 18, 1161–1168 (2021). https://doi.org/10.1038/s41592-021-01254-9
