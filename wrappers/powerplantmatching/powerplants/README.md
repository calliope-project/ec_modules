# Wrapper for powerplantmatching powerplants

Runs `powerplantmatching`'s matching process with optional custom configurations.

```mermaid
flowchart LR
    I1(config.yaml) --> |Optional|W((powerplants))
    W --> O1(powerplants.csv)
    W --> |Optional| O2(stats_aggregated.png)
```

>[!warning]
>Although useful, `powerplantmatching` is a risky datasource for the following reasons:
>
>- Updated data tends to change between releases due to the unstable nature of the algorithm employed.
>- DOUBLE COUNTING MAY STILL HAPPEN!
>- The resulting data is not a perfect match with national statistics.
>
>Use with care!
