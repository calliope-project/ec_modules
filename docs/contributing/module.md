# Module requirements

## File structure

We attempt to follow `snakemake`'s [recomended structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility) to ensure our modules can be easily exported to other projects.
Additionally, we extend it with a few quality of life improvements.

??? example "Example of a fully featured module"
    Most modules should follow a structure similar to the one below:
    ```tree
    modules/example/
    ├── AUTHORS
    ├── LICENSE
    ├── README.md
    ├── rulegraph.png
    ├── config/
    │   └── default.yaml
    ├── docs/
    │   └── extra_docs.md
    ├── resources/
    │   ├── download1.nc
    │   └── download2.csv
    ├── results/
    │   ├── plots/
    │   │   └── graph.png
    │   ├── table_data.csv
    │   └── spatial_data.nc
    └── workflow/
        ├── Snakefile
        ├── envs/
        │   ├── geo.yaml
        │   └── shell.yaml
        ├── internal/
        │   ├── settings.yaml
        │   └── mapping.csv
        ├── profiles/
        │   └── default/
        │       └── config.yaml
        ├── report/
        │   └── review.rst
        ├── rules/
        │   ├── rules1.smk
        │   └── rules2.smk
        ├── schemas/
        │   └── default.schema.yaml
        └── scripts/
            ├── example.py
            └── README.md
    ```

### Obligatory components

Please ensure that your module has the following:

- An `AUTHORS` file. This lists all the people who have contributed to a module.
- A `LICENSE` file. This should refer to the `AUTHORS` file above. We recommend to use the [MIT license](https://opensource.org/license/mit).
- A `README.md` file. Here you describe the module's functionality (see [documentation](#documentation)).
- A `config/default.yaml` file. This contains the module's default configuration.
It serves only as an example: users will override it using `snakemake`'s [`module` functionality](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
- A `resources/` folder.  All **downloaded data** used by your workflow should be deposited here.
- A `results/` folder. All **rule outputs** of your workflow should be deposited here.

    ??? info
        There is no real difference between `resources/` and `results/` in terms of `snakemake` functionality.
        This separation is mostly to enable easier bookkeeping and to help module users differentiate between the data you _downloaded_ and data you have _processed_.

- A `workflow/` folder with the following:
    - A `Snakefile`. This contains the module's main `snakemake` functionality (see [ensuring modularity](#ensuring-modularity)).
    - An `envs/` folder. All your `conda` environments live here.
    - An `internal/` folder. Small files needed by your workflow should be placed here.
    These can be `.yaml` files with internal configuration or `.csv` files with mappings to aid in parsing.
    - A `profiles/default/config.yaml` file. Our standard `snakemake` [profile](https://snakemake.readthedocs.io/en/v8.18.0/executing/cli.html#profiles) configuration.

        ??? warning
            `profiles/default/config.yaml` should not be altered!
            This file's settings ensure that wrappers and conda environments execute seamlessly across modules.

    - A `rules/` folder. This contains all the `rule.smk` files used by your module.
    - A `scripts/` folder. All python / R code used by your rules should be here.

### Optional components

You can additionally compliment your module with the following:

- A `docs/` folder. If needed, additional documentation can be placed here.
- In `workflow/`:
    - A `reports/` folder. Automatic `snakemake` [reports](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html) can be placed here. This will allow users to easily follow your workflow, evaluate runtimes and even visualise plots!
    - A `schemas/` folder. Schemas used [validate](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#validation) user configuration should be placed here.

## Documentation

The `README.md` file should be a pragmatic quick example of what your module needs to function.
It should contain at least the following things:

- [X] A mermaid diagram briefly showing its inputs and outputs.
- [X] A DAG showing the rule order.
- [X] A citation text.
- [X] A references section (if your workflow is based on the work of others).

If additional context is needed, please place it in the `docs/` folder.

??? example "Example: Hydropower module"
    **Euro-Calliope hydropower data**

    Easily generate timeseries for hydropower plants for any region in Europe.

    **_Input-output_**

    ```mermaid
    flowchart LR
        id1[("
            ERA5
            HydroBASINS
        ")] --> |Download| M
        C[/"
            config.yaml
            - shapefile.geojson
            - powerplants.csv
        "/] -->M((hydropower))
        M --> O1("
            capacity-factors-RoR.csv
            capacity-factors-basins.csv
            ")
        M --> O2("
            region-power-capacity.csv
            region-storage-capacity.csv
            ")
    ```

    **_DAG_**

    ![DAG](images/rulegraph.png)

    **_Citation_**

    Tröndle, T., & Pickering, B. (2021). Euro-Calliope Hydropower data [Computer software]. https://doi.org/10.5281/zenodo.3949793

    **_References_**

    None

## Ensuring modularity

Modules are essentially workflows that can be exported to other projects. They operate largelly the same way, with some key differences:

- Plain string references to files _outside_ of `rules` will not work. This happens because `snakemake` will use the path of the importing workflow as the base path. To ensure this does not cause issues, always use `workflow.source_path()` when requesting files in folders like `internal/`, `schema/` etc.

    ??? example "Example: loading a schema"
        :cross_mark: This will fail...

        ```python
        validate(config, "schema/config.schema.yaml")
        ```
        :white_check_mark: And this will work!

        ```python
        validate(config, workflow.source_path("schema/config.schema.yaml"))
        ```

- References to files _within_ rules will work the same way they do in regular workflows as long as you avoid direct paths. No need to use special functions here!

    ??? example "Example: defining a rule"
        :white_check_mark: This will work just fine!
        Even `shadow`, which uses a temporary directory.

        ```smk
        rule raw_population_unzipped:
            message: "Extract the population data TIF."
            input: rules.download_raw_population_zipped.output
            shadow: "minimal"
            output: "results/population.tif"
            conda: "../envs/shell.yaml"
            shell: "unzip -p '{input}' 'JRC_1K_POP_2018.tif' > '{output}'"
        ```

- The configuration in `config/default.yaml` will _always_ be overriden by users. This makes documentation quite important, as `default.yaml` configuration is just a suggestion!
- Both module developers and users must use our standard profile in `workflow/profiles/default/config.yaml`. Otherwise, wrappers and `conda` won't work as intended.
