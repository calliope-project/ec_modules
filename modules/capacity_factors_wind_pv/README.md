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
        capacityfactors (.nc)
        eez (.gpkg)
        ...
    `")] --> |Download| M
    C1[/"`**User input**
        units.geojson
        ...
    `"/] --> |Resources| M((capacity_factors_wind_pv))
    M --> O1("
        capacityfactors-wind-offshore.csv
        ")
    M --> O2("
        capacityfactors-wind-onshore.csv
        ")
    M --> O3("
        capacityfactors-wind-offshore.csv
    ")
    M --> O4("
        capacityfactors-open-field-pv.csv
    ")
    M --> O5("
        capacityfactors-rooftop-pv.csv
    ")
    M --> O6("
        capacityfactors-rooftop-pv-n.csv
    ")
    M --> O7("
        capacityfactors-rooftop-pv-e-w.csv
    ")
    M --> O8("
        capacityfactors-rooftop-pv-s-flat.csv
    ")
```

## DAG

Here is a brief example of the module's steps.

![DAG](rulegraph.png)
