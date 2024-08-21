"""Capacity layout wrapper."""

import atlite
import xarray as xr
from matplotlib import pyplot as plt
from pyproj import CRS
from shapely import box

M2_TO_KM2 = 1 / 1e6

availability = xr.open_dataarray(snakemake.input.availability_matrix)
cutout = atlite.Cutout(
    path=snakemake.input.cutout, **snakemake.params.get("cutout_kwargs", {})
)
assert CRS.from_user_input(cutout.crs).is_geographic
assert box(*cutout.bounds).covers(
    box(
        float(availability.x.min()),
        float(availability.y.min()),
        float(availability.x.max()),
        float(availability.y.max()),
    )
)

layout_crs = snakemake.params.layout_projected_crs
assert CRS.from_user_input(layout_crs).is_projected
layout = (
    cutout.uniform_density_layout(snakemake.params.max_units_per_km2, crs=layout_crs)
    * M2_TO_KM2
    * availability.sum(set(availability.dims) - {"x", "y"})
)
layout.name = snakemake.params.layout_name
layout.attrs["units"] = snakemake.params.units
layout.to_netcdf(snakemake.output.layout)

plot_layout = snakemake.output.get("plot_layout")
if plot_layout:
    fig, ax = plt.subplots()
    layout.plot()
    cutout.grid.plot(ax=ax, color="None", edgecolor="grey", ls=":")
    plt.savefig(plot_layout)
