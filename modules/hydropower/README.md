# Module: Hydropower

A module to estimate hydropower capacity factors and powerplant capacities for arbitrary regions in Europe.

## Input-Output

Here is a brief summary of the IO structure of the module.

Users must specify a remote location of shapefile with the desired subregions in the configuration, which will be downloaded and processed into timeseries and capacity values by the module.

```mermaid
flowchart LR
    I1(shapefile.geojson) -.-> |Download| C
    C(config.yaml) -->M((hydropower))
    M --> O1(capacity-factors-RoR.csv)
    M --> O2(capacity-factors-basin.csv)
    M --> O3(region-power-capacity.csv)
    M --> O4(region-storage-capacity.csv)
```

## DAG

Here is a brief overview of the module's steps.
Please consult the code for more details.

![dag](rulegraph.png)
