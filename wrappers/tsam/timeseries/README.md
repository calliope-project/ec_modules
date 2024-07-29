# Wrapper for `tsam` time series aggregation methods

Runs `tsam` time series aggregation methods for a list of time series files.
Optionally, outputs accuracy indicator data and index matching data.

>[!warning]
>It is assumed that all time series is in the same datetime format and has matching rows!

## Example

```snakemake
rule tsam_timeseries:
    input:
        "timeseries1.csv",
        "timeseries2.csv"
    output:
        typical_periods = "output/typical_periods.csv",
        accuracy_indicators = "output/accuracy_indicators.csv",
        index_matching = "output/index_matching.csv"
    params:
        noTypicalPeriods = 8,
        hoursPerPeriod = 24,
        segmentation = True,
        noSegments = 8,
        representationMethod = "distributionAndMinMaxRepresentation",
        distributionPeriodWise = False,
        clusterMethod = 'hierarchical'
    threads: 4
    wrapper: github("calliope-project/ec_modules", path="wrappers/tsam/timeseries")
```
