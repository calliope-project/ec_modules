# Easy Energy Modules - demand_electricity

A module preparing electricity demand time series

## Input-Ouput

Here is a brief IO diagram of the module's operation.

```mermaid
---
title: demand_electricity
---
flowchart LR
    D1[("`**Databases**
        raw-potentials/demand.csv
        raw-load-data.csv
        ...
    `")] --> |Download| M
    C1[/"`**User input**
        units.geojson
        ...
    `"/] --> |Resources| M((demand_electricity))
    M --> O1("
        demand-electricity-ehighways.csv
        ")
    M --> O2("
        electricity-demand-national.csv
        ")
```

## DAG

Here is a brief example of the module's steps.

![DAG](rulegraph.png)
