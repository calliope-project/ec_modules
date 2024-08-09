# `atlite` land-use availability

For the given shapes and cutout, combine raster files to filter land area. The result is a netCDF file with the usable land ratio per cell.

>[!important]
>`raster`, `raster_codes` and `raster_kwargs` must always match in length!
>Use `{}` if no kwargs want to be passed.

```mermaid
flowchart LR
    I1(cutout.nc) --> W((availability))
    I2(shapefile.shp) --> W
    I3 --> W
    I3 --> W
    I3(raster1.tif, raster2.tif...) --> W
    W --> O1(availability_matrix.nc)
    W --> |Optional| O2(plot_shape_availability.png)
    W --> |Optional| O3(plot_availability_matrix.png)
```
