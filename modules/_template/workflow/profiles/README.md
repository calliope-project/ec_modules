# Profiles

>[!WARNING]
>Do not modify the files in this folder

`default/config.yaml` is a `snakemake` [profile](https://snakemake.readthedocs.io/en/v8.18.0/executing/cli.html) that enables us to link `ec_modules` together.
Any workflow that wishes to use our modules and wrappers should have a similar file.
It ensures two things:

- That the workflow is always ran using `conda`.
- That wrappers point to the [`ec_modules`](https://github.com/calliope-project/ec_modules) repository.
