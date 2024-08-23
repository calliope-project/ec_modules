# The easy energy system workflows and wrappers repository

This aims to be a collection of helper workflows and wrappers that aid in quickly developing Energy System (ESM) models, inspired by Snakemake's [workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/) and [wrapper repository](https://snakemake-wrappers.readthedocs.io/en/stable/).
Its main purpose is enabling researchers to share data workflows between studies to avoid the bloat problem in large energy system workflows.

## Using our modules

Please follow our user guidelines to get quickly set-up to re-using our modelling tools!

## Development setup

If you wish to test the whole repository, please run.

```shell
mamba env create -f environment.yaml
mamba activate ec_modules
```

>[!warning]
> The `ec_modules` environment is mostly used for development and testing!
