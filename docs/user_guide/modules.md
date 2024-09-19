# Using Modules

Think of them as **workflows that can be exported to other projects**.
Their settings can be **re-configured**, allowing you to reprocess data with different parametric assumptions.
Some use cases:

- Aggregating transmission line data into regions.
- Calculating renewable potentials for a given set of regions.
- Combining datasets of existing power production facilities.
- And many more!

Generally, modules are used in the following steps:

1. You configure the module to behave a certain way.
2. You give the module some files to process.
3. You collect the results.

The module will automatically download necessary files and execute code to produce the requested output.

You can read more about modules in the [`snakemake` documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).

## What is a module?

Modules can be visualised as a DAG (direct acyclic graph) with a series of steps that "connects" to your workflow.

??? example "Module as a DAG"

    The module (blue) connects to your general workflow (yellow).

    ![Module example](images/module.png)

However, it's generally more useful to think of modules as an IO (input-output) diagram.
After they are given certain inputs and configuration, they will automatically process them and produce output files.

??? example "Module as IO"

    In this case, a hydropower module needs two **user** inputs (`regions.geojson`, `powerplants.csv`) and has four outputs.

    ```mermaid
    flowchart LR
        id1[("`**Automatic**
            ERA5
            HydroBASINS
        `")] --> |Download| M
        C[/"`**User**
            regions.geojson
            powerplants.csv
        `"/] -->M((hydropower))
        M --> O1("`**Timeseries**
            capacity-factors-RoR.csv
            capacity-factors-basins.csv
            `")
        M --> O2("`**Installed capacity**
            region-power-capacity.csv
            region-storage-capacity.csv
            `")
    ```

## Configuring a module

You can always find a default configuration in the following location:

```tree
module_name/
├── README.md
└── config/
    └── default.yaml      # Here!
```

The module's schema files (`module/workflow/schemas/`) will have detailed information on how to configure the module.
The number and purpose of the configuration options will vary per module.

We recommend storing module configuration files separately, and adding them to your general configuration with module-specific keys.
This will ensure there are no conflicts between modules and your general configuration.

??? example "Module configuration example:"

    This is how you should arrange your files:

    ```tree
    your_workflow/
    ├── README.md
    └── config/
        ├── general.yaml            # Your general configuration!
        └── modules/
            ├── heat.yaml           # Module-specific configuration!
            ├── hydropower.yaml
            └── transport_road.yaml
    ```

    And this is how a `hydropower.yaml` above would look like:

    ```yaml
    module-hydro:
        use_default_user_resources: False
        hydro_power_database:
            version: "v10"
        HydroBASINS_database:
            level: "07"
    ```

## Using our modules

### Import

Adding a module to your `snakemake` workflow is easy, as long as you have followed our [getting started](getting_started.md) instructions.
All you need to do is import it with a standardised structure.

???+ example "Using a `hydropower` module"

    In this example, assume you wish to obtain hydropower data.
    You must tell `snakemake` the location of the module in GitHub, and pass your configuration to it.

    ```python
    # Extend your general configuration with your module-specific configuration
    configfile: "config/modules/hydropower.yaml"

    module hydropower:
        # Request snakemake to use a specific module version (tag)
        snakefile:
            github(
            "calliope-project/ec_modules",
            path="modules/hydropower/Snakefile",
            tag="v1.0.0"
            )
        # Provide your module-specific configuration
        config: config["hydropower"]
        # Add a prefix to avoid file conflicts
        # (e.g., results/output.csv -> module_hydropower/results/output.csv)
        prefix: "module_hydropower"

    # Rewrite rule names to avoid conflicts (e.g., all -> module_hydro_all)
    use rule * from hydropower as module_hydropower_*
    ```

??? example "using `prefix:`"

    `prefix:` enables you to avoid file name conflicts between modules by prepending all paths in the module with an additional folder.
    This also affects the paths in your configuration!
    In the [configuration example](#configuring-a-module) above, `"resources/shapes_spain.geojson"` will _also_ be affected!

    So, you should place your `shapes_spain.geojson` file in the following location:

    `prefix:` + `resources/shapes_spain.geojson` -> `module_hydropower/resources/shapes_spain.geojson`

??? warning "`snakemake` compatibility"

    Some modules might be incompatible with older versions of `snakemake`.
    We recommend to update `snakemake` to the newest version possible in such cases.
    We generally enforce `snakemake = 8.10` as the minimum supported version.

### Adding input files

Afterwards, all you need to do is place input files in the module's `resources/user` folder (either manually or via `snakemake`).
Each module will specify the file it needs in its documentation, and provide examples.

??? example "Inputting a national resolution file"

    Our example module needs a file in the format `shapes_{resolution}.geojson`, which can contain an arbitrary regional resolution.

    Assuming that you want data aggregated at a national level, you should use a `.geojson` file with that resolution as input.

    ```python
    rule hydropower_input:
        message: "Input the desired resolution to the module."
        input: "data/national_resolution.geojson"
        output: f"module_hydropower/resources/user/shapes_national.geojson"
        conda: "../envs/shell.yaml"
        shell: "cp {input} {output}"
    ```

### Requesting results

To obtain results, all you need to do is tell `snakemake` to output a specific file that the module produces.
Examples for all relevant files will be detailed in the module's documentation.

??? example "Requesting national data"

    After placing your spatial resolution file, all you need to do is request the timeseries file, respecting the same wildcards.

    In this case, the module specifies

    `results/{resolution}/{year}/run_of_river_timeseries.csv`

    Assuming you want data for 2020, you would do the following:

    ```python
    rule hydropower_output:
        message: "Request hydropower timeseries."
        # Request timeseries for 2020 specifically.
        input: f"module_hydropower/results/national/2020/run_of_river_timeseries.csv"
        # Copy it to a different location for further processing.
        output: f"build/model/timeseries/demand/electricity.csv"
        conda: "../envs/shell.yaml"
        shell: "cp '{input}' '{output}'"
    ```
