# Easy Energy Modules - biofuels

A module preparing biofuels cost and potentials

## Input-Ouput

Here is a brief IO diagram of the module's operation.

```mermaid
---
title: biofuels
---
flowchart LR
    D1[("`**automatic**
        biomass potentials (JRC ENSPRESO)
        population (JRC GHS-POP)
        land cover (https://doi.org/10.5281/zenodo.3244985)
        ...
    `")] --> |Download| M
    C1[/"`**user**
        spatial_units.geojson
        ...
    `"/] --> M((biofuels))
    M --> O1("
        results/{resolution}/{scenario}/potential-mwh-per-year.csv
        ")
    M --> O2("
        results/{resolution}/{scenario}/costs-eur-per-mwh.csv
        ")
```

User
----
- resources/user/spatial_units.geojson

Output
------
- results/ehighways/high/potential-mwh-per-year.csv
- results/ehighways/high/costs-eur-per-mwh.csv
- results/ehighways/low/potential-mwh-per-year.csv
- results/ehighways/low/costs-eur-per-mwh.csv
- results/ehighways/medium/potential-mwh-per-year.csv
- results/ehighways/medium/costs-eur-per-mwh.csv

## DAG

Here is a brief example of the module's steps.

![DAG](rulegraph.png)
