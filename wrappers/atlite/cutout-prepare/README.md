# `atlite` cutout-prepare

Download cutout data for the requested shapefile, prepare it for `atlite` functions, and save it locally for reuse.

>[!important]
>`atlite` may ignore parameters if the cutout already exists!

```mermaid
flowchart LR
    I1(shapefile.shp) --> W((cutout-prepare))
    I2(gebco.nc) --> |Optional| W
    W --> O1(cutout.nc)
    W --> |Optional|O2(plot_cutout.png)
```
