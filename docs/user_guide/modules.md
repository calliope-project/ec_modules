# [Modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules)

Think of them as **workflows that can be exported to other projects**.
Their settings can be **re-configured**, allowing you to reproduce data with different parametric assumptions.
They have one disadvantage: their **inputs are static**, meaning that you cannot pass results from other snakemake rules.
Modules are best for stable workflows whose inputs are not expected to change.
Some use cases:

- Converting transmission line data into nodes / regions.
- Calculating renewable potentials for a given set of regions.
- Combining datasets of existing power production facilities.
- And many more!

## Visualising a module

Modules can be visualised as a regular DAG (direct acyclic graph) that "connects" to your workflow.

??? example "Module as DAG"
    ![Module example](images/module.png)

However, it's generally more useful to think of them as an IO (input-output) diagram that you can influence by reconfiguring it.
The databases used and the order of the processing steps may not change, but the output can be influenced depending on the configuration files you give it.

!!! example "Module as IO"

    In this case, the module's results will change depending on the shapes and powerplants that the configuration points to.

    ??? info "Hydropower workflow"

        ![Hydropower workflow](images/module_rulegraph.png)

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

## Importing a module

Adding a module to your `snakemake` workflow is as easy as:

```python
module hydropower:
    # Plain paths, URLs and github / gitlab calls are possible
    snakefile:
        github(
          "calliope-project/ec_modules",
          path="modules/hydropower/Snakefile",
          tag="v1.0.0"
        )
    # Provide your own configuration to the module
    config: config["hydropower"]
    # Add a prefix to the results to avoid file conflicts
    # (e.g., results/output.csv -> module-hydropower/results/output.csv)
    prefix: "module-hydropower"

# Rewrites rule names to avoid conflicts (e.g., all -> module_hydro_all)
use rule * from hydropower as module_hydropower_*
```

## Configuring a module


## Using a module
