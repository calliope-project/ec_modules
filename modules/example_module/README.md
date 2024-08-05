# Module example

General file structure to follow for `snakemake` modules.

>[!important]
>Always call snakemake at the `example_module/` level, not at the workflow level!
>This is by `snakemake` design.

```ascii
example_module/
┣ config/                # Default configuration, can be overridden!
┃ ┗ config.yaml
┣ resources/             # Static files needed by your workflow
┣ results/               # Put all rule outputs here!
┣ workflow/
┃ ┣ envs/                # Conda environments
┃ ┃ ┗ example_env.yaml
┃ ┣ report/              # For snakemake's report functionality
┃ ┣ rules/               # Rule files
┃ ┃ ┗ example.smk
┃ ┣ schemas/             # Schemas to check configuration files
┃ ┃ ┗ config.schema.yaml
┃ ┣ scripts/             # Actual code!
┃ ┃ ┗ example.py
┃ ┗ Snakefile            # main rule lives here!
┣ LICENSE
┗ README.md
```
