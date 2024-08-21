# `atlite` pv

Uses a cutout with pv-relevant features to produce a dataset with:

- Profile time-series (either capacity factor or maximum generation)
- Maximum installable-capacity per region

Optionally, get plots for the generation profile and max. capacity.

```mermaid
flowchart LR
    I1(cutout.nc) -->W((pv))
    I2(shapes.shp) --> W
    I3(layout.nc) --> W
    W --> O1(dataset.nc)
    W -->|Optional|O2(plot_profile.png)
    W -->|Optional|O3(plot_max_capacity.png)
```
