# Easy Energy Modules - industry

A module to estimate industrial energy demand and decarbonisation potential in European subregions.

## Input-Ouput

Here is a brief IO diagram of the module's operation.

```mermaid
---
title: industry
---
flowchart LR
    D1[("`**Databases**
        eurostat
        Swiss BFE
        JRC-IDEES
    `")] --> |Download| M
    C1[/"`**User input**
        shapefile.geojson
    `"/] --> |Resources| M((industry))
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
