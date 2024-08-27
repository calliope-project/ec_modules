# Easy Energy Modules - capacity_factors_wind_pv

A module that prepares capacity factors for pv and wind, both onshore and offshore

## Input-Ouput

Here is a brief IO diagram of the module's operation.

```mermaid
---
title: capacity_factors_wind_pv
---
flowchart LR
    D1[("`**Databases**
        database1
        ...
    `")] --> |Download| M
    C1[/"`**User input**
        shapefile.geojson
        ...
    `"/] --> |Resources| M((capacity_factors_wind_pv))
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
