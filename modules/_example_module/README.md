# Module example

Please follow this file structure when developing modules.
It is based on `snakemake`'s [recommended structure to ensure reproducibility](https://snakemake.readthedocs.io/en/v8.18.0/snakefiles/deployment.html), and will save you a lot of trouble.

>[!important]
>Always call snakemake at the `example_module/` level:
>
>`ec_modules/modules_example_module$ snakemake -c 1`

```text
_example_module/
├── config
│   └── default.yaml             # Default configuration. Can be overridden by users.
├── LICENSE
├── results                      # Place all your output files here.
└── workflow
    ├── envs
    │   └── shell.yaml           # All rules should use a conda environment defined here.
    ├── profiles
    │   └── default
    │       └── config.yaml      # Do not delete!
    ├── report                   # Optional, generate automatic reports with snakemake!
    ├── resources
    │   └── internal.yaml        # Small files needed by your module, such as internal configuration.
    ├── rules
    │   └── example.smk          # Place all your workflow rules here.
    ├── schemas
    │   └── config.schema.yaml   # Optional, validate user configurations automatically!
    ├── scripts
    │   └── example.py           # All your scripts should go here.
    └── Snakefile
```
