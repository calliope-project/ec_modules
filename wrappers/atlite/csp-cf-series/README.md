# Concentrated Solar Power capacity factor series with `atlite`

Uses a cutout with CSP relevant data to produce capacity factor timeseries for energy system models.
Optionally, get some nice plots with the average CF per calculated region.

```mermaid
flowchart LR
    I1(Cutout.nc) -->W((csp-cf-series))
    I2(Shapefile.shp) --> W
    W --> O1(timeseries.csv)
    W -->|Optional|O2(plot_mean_cf.png)
```

## Example

```snakemake
code here
```
