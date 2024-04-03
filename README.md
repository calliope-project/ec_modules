# Snakemake modularisation tests

This toy repo contains trials of how to modularise Snakemake repositories.

Currently, [Sector-Coupled Euro Calliope (SCEC)](https://github.com/calliope-project/sector-coupled-euro-calliope) and [Euro-Calliope (EC)](https://github.com/calliope-project/euro-calliope) are quite user unfriendly in their Snakemake steps. EC has more than 60 rules, and SCEC adds more than 40 on top of it!

**Modularisation is needed to abstract complexity, letting people focus on doing science instead!**

## Setup

Please run the following commands:

```shell
mamba env create -f environment.yaml
mamba activate ec_modules
```

## Goals

Successful modularisation must meet the following criteria to ensure maximum compatibility and keep conflicts at a minimum.

- **Inputs must be file/folder paths**: to keep modules agnostic to project structure.
- **Outputs should only be files**: meaning that code variables cannot be passed between modules.
- **Outputs must be specified by templates**: modules might need to use the output of other modules as input (e.g., a _pv_capacity_ module might use a `*.geojson` file produced by a _regions_ module). To avoid breaking functionality, output template files must be specified in the `templates/` folder of a module to let other developers know what they'll be working with.
- **File structure should be standardised**: to help users reach relevant files quickly, modules must at least have _src_, _out_, _templates_ and _tmp_ folders.
- **Containerised**: each module must use its own _env_ to keep version conflicts at a minimum.
- **Isolated configuration YAML**: some sort of module configurability is needed to keep things flexible. However, modules must only have access to their own configuration parameters to avoid potential conflicts. Providing template configuration files is also recommended.
- **Debuggable and testable**: problems are bound to happen, so modules should be relatively easy to debug, and should include a few test cases by default.

<p align="center">
  <img src="docs/images/module_setup.png" />
</p>

A detailed summary of potential modularisation candidates can be found [here](docs/euro-calliope%20DAG%20structure%202024-04-02.pdf).

## Proposals

These are tryouts since Snakemake is very picky on how code is included.

To do:

- [x] Find a pseudo-modularisation method within a single repo (most friendly to current calliope approaches).
- [ ] Find an approach using Snakemake's wrapper functionality.
- [ ] If not all requirements are fulfilled, look for alternative approaches.
