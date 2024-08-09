"""atlite technology capacity factor and maximum capacity wrapper.

Only valid for certain technologies.
"""

import atlite
import geopandas as gpd
import xarray as xr
from matplotlib import pyplot as plt
from pyproj import CRS
from shapely import box

VALID_TECHS = ("pv", "wind", "csp")

cutout = atlite.Cutout(
    path=snakemake.input.cutout, **snakemake.params.get("cutout_kwargs", {})
)
assert CRS.from_user_input(cutout.crs).is_geographic

tech = snakemake.params.tech
assert tech in VALID_TECHS
tech_kwargs = snakemake.params.tech_kwargs

shapes = gpd.read_file(snakemake.input.shapefile)
assert shapes.crs.is_geographic
assert box(*cutout.bounds).covers(box(*shapes.total_bounds))
shapes = shapes.set_index(snakemake.params.shapefile_name_column)
tech_kwargs["shapes"] = shapes

layout = xr.open_dataarray(snakemake.input.layout)
assert box(*cutout.bounds).covers(
    box(
        float(layout.x.min()),
        float(layout.y.min()),
        float(layout.x.max()),
        float(layout.y.max()),
    )
)
tech_kwargs["layout"] = layout

tech_function = getattr(cutout, tech)
profile, capacity = tech_function(return_capacity=True, **tech_kwargs)

data = xr.Dataset()
data["profile"] = profile
data["max_nominal_capacity"] = capacity
data.to_netcdf(snakemake.output.dataset)
if snakemake.output.get("plot_profile"):
    fig, ax = plt.subplots()
    profile.to_pandas().plot(legend=False, title=tech)
    plt.tight_layout()
    plt.savefig(snakemake.output.plot_profile)
if snakemake.output.get("plot_max_capacity"):
    fig, ax = plt.subplots()
    capacity.to_pandas().plot.bar(rot=45, ylabel=capacity.attrs["units"], title=tech)
    plt.tight_layout()
    plt.savefig(snakemake.output.plot_max_capacity)
