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
rule atlite_csp_cf_series:
    input:
        cutout = "cutout_csp.nc",
        shapefile = "portugal.geojson"
    output:
        timeseries = "output/portugal.csv",
        plot_mean_cf = "output/mean_cf.png"
    params:
        shapefile_name_column = "state",
        installation = "SAM_solar_tower",
    threads: 4
    wrapper: github("calliope-project/ec_modules", path="wrappers/atlite/csp-cf-series")
```
