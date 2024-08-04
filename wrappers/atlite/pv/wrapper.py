"""atlite PV capacity factor wrapper."""

import atlite
import geopandas as gpd
import xarray as xr
from matplotlib import pyplot as plt
from pyproj import CRS
from shapely import box

cutout = atlite.Cutout(
    path=snakemake.input.cutout, **snakemake.params.get("cutout_kwargs", {})
)
assert CRS.from_user_input(cutout.crs).is_geographic

pv_kwargs = snakemake.params.pv_kwargs

shapes = gpd.read_file(snakemake.input.shapefile)
assert shapes.crs.is_geographic
assert box(*cutout.bounds).covers(box(*shapes.total_bounds))
shapes = shapes.set_index(snakemake.params.shapefile_name_column)
pv_kwargs["shapes"] = shapes

layout = xr.open_dataarray(snakemake.input.layout)
assert box(*cutout.bounds).covers(
    box(
        float(layout.x.min()),
        float(layout.y.min()),
        float(layout.x.max()),
        float(layout.y.max()),
    )
)
pv_kwargs["layout"] = layout

profile, capacity = cutout.pv(return_capacity=True, **pv_kwargs)
data = xr.Dataset()
data["profile"] = profile
data["max_nominal_capacity"] = capacity
data.to_netcdf(snakemake.output.dataset)
if snakemake.output.get("plot_profile"):
    fig, ax = plt.subplots()
    profile.to_pandas().plot(legend=False, title="pv")
    plt.savefig(snakemake.output.plot_profile)
if snakemake.output.get("plot_max_capacity"):
    fig, ax = plt.subplots()
    capacity.to_pandas().plot.bar(rot=45, ylabel=capacity.attrs["units"], title="pv")
    plt.savefig(snakemake.output.plot_max_capacity)
