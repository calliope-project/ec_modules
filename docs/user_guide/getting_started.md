# Getting started

In this section we aim to get you familiarised with the `ec_modules` approach as quickly as possible.
More in-dept information is available the `snakemake` documentation, but it's generally [quite](https://snakemake-wrappers.readthedocs.io/en/stable/) [spread](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html) [out](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).
To make things easier, we summarize relevant functionality here.

## Pre-requisites

In order to run our workflows, you need `conda`/`mamba` environment with `snakemake`.
Feel free to skip this if you already have a workflow and just wish to integrate our modules into it.

1. Install `mamba` or `conda`. We recommend following `mamba`'s [installation advice](https://github.com/mamba-org/mamba), as it is much faster than regular `conda`.
2. Create an environment with `snakemake>=8.10` installed.

    `mamba env create -n my_workflow bioconda::snakemake-minimal>=8.10`

3. Create a workflow using `snakemake`'s [standard template](https://github.com/snakemake-workflows/snakemake-workflow-template).

## Configuring your workflows to use `ec_modules`

Adding compatibility with `ec_modules` is simple!
Just follow these two steps.

1. Create a [default configuration profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) in `your_project/workflow/profiles/default/config.yaml`.

    ??? example
        For `your_project/`:

        ```yaml
        your_project/
        ├── environment.yaml
        ├── LICENSE.md
        ├── README.md
        ├── config/
        └── workflow/
            ├── Snakefile
            ├── config/
            ├── envs/
            └── profiles/
                └── default/
                    └── config.yaml   # Here!
        ```

2. Add the following configuration:

    ```yaml
    software-deployment-method: conda
    use-conda: True
    wrapper-prefix: https://github.com/calliope-project/ec_modules/raw/
    ```

That's it!
You're ready to use reproducible energy model workflows.
