# `tsam` - time series aggregation methods

Runs `tsam` time series aggregation methods for a list of time series files.
Optionally, outputs accuracy indicator data and index matching data.

>[!INFO]
>Multiple variables per file are possible, but the variable names must be unique across all files!

```mermaid
flowchart LR
    I1(data1.csv) --> W((timeseries))
    I2(data2.csv) --> W
    I3(data3.csv) --> W
    W --> O1(typtical_periods.csv)
    W --> |Optional| O2(accuracy_indicators.csv)
    W --> |Optional| O3(index_matching.csv)
```
