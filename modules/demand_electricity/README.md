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
        database1
        ...
    `")] --> |Download| M
    C1[/"`**User input**
        shapefile.geojson
        ...
    `"/] --> |Resources| M((demand_electricity))
    M --> O1("
        output1.csv
        ")
    M --> O2("
        output2.nc
        ")
```

## DAG

Here is a brief example of the module's steps.

![DAG](rulegraph.png)
