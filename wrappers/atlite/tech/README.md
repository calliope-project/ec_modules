# `atlite` technology wrapper

Uses a cutout with technology-relevant features to produce a dataset with:

- Profile time-series (either capacity factor or maximum generation)
- Maximum installable-capacity per region

Optionally, get some nice plots for the profile and max. capacity.

>[!important]
>This helper wrapper is only valid for the following technology functions in `atlite`: `pv`, `wind`, `csp`.
>
>Keep in mind that the required parameters change between them!

```mermaid
flowchart LR
    I1(cutout.nc) -->W((tech))
    I2(shapes.shp) --> W
    I3(layout.nc) --> W
    W --> O1(dataset.nc)
    W -->|Optional|O2(plot_profile.png)
    W -->|Optional|O3(plot_max_capacity.png)
```
