# Wrapper for `tsam` time series aggregation methods

Runs `tsam` time series aggregation methods for a list of time series files.
Optionally, outputs accuracy indicator data and index matching data.

>[!warning]
>Please ensure that all timeseries have the following format:
>|                     | T    | Load        |
>|---------------------|------|-------------|
>| 2009-12-31 23:30:00 | -2.1 | 375.4783938 |
>| 2010-01-01 00:30:00 | -2.8 | 364.5413263 |
>| 2010-01-01 01:30:00 | -3.3 | 357.4168443 |
>| 2010-01-01 02:30:00 | -3.2 | 350.1913058 |
>
>Simiarly, multiple variables per file are possible. However, the same variable names must be unique across all files!

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
