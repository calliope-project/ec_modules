<!-- Please provide a concise summary of the module in this section. -->
<!-- --8<-- [start:intro] -->
# demand_electricity

A module preparing electricity demand time series for European nations.

<!-- --8<-- [end:intro] -->

## Input-Ouput

<!-- Please fill in this diagram including: wildcards, user resources and final results. -->
<!-- --8<-- [start:mermaid] -->

```mermaid
---
title: demand_electricity
---
flowchart LR
    D1[("`**automatic**
        raw-potentials/demand.csv
        OPSD - time_series
        ...
    `")] --> |Download| M
    C1[/"`**eser**
        shapes_{resolution}.geojson
        ...
    `"/] --> M((demand_electricity))
    M --> O1("`**timeseries**
        {resolution}/{year}/demand_electricity.csv
        `")
    M --> O2("`**plots**
        {resolution}/{year}/plot_map.png
        {resolution}/{year}/plot_timeseries.png
        `")
```
<!-- --8<-- [end:mermaid] -->

### Wildcards

<!-- Please explain what wildcards are required by users here. -->
<!-- --8<-- [start:wildcards] -->
- **{resolution}**: Determines the number of regions that the module will process. Importantly, it must be specified for the correct input file, which can be obtained from [Euro-Calliope datasets](https://zenodo.org/records/6600619). The following options are possible:
    - national
    - regional
    - ehighways
    - continental
- **{year}**: The year for which historical demand will be generated. Valid values: 2005 up to 2019.
<!-- --8<-- [end:wildcards] -->

### User

<!-- Please briefly explain user resources here. -->
<!-- --8<-- [start:user] -->

- **resources/user/shapes_{resolution}.geojson**: a file with the shapes in the desired spatial resolution. Can be obtained from [Euro-Calliope datasets](https://zenodo.org/records/6600619).

<!-- --8<-- [end:user] -->

### Results

<!-- Please briefly explain final result files here. -->
<!-- --8<-- [start:results] -->

- **{resolution}/{year}/demand_electricity.csv**: annual electricity demand timeseries per region.
- **{resolution}/{year}/plot_map.png**: map showing the distribution of annual demand.
- **{resolution}/{year}/plot_timeseries.png**: annual demand timeseries and load duration curve.

<!-- --8<-- [end:results]  -->

## References

<!-- Please cite studies and datasets used for this workflow below. -->
<!-- --8<-- [start:references] -->

- Open Power System Data (2020). Time series data package [Dataset]. <https://github.com/Open-Power-System-Data/time_series>

<!-- --8<-- [start:references] -->
