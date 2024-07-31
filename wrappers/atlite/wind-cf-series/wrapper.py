"""Wind capacity factor wrapper for atlite."""

import atlite
import geopandas as gpd
from matplotlib import pyplot as plt

shapes = gpd.read_file(snakemake.input.shapefile)

column = snakemake.params.shapefile_name_column
shape_names = snakemake.params.get("shapefile_specific_names")
if shape_names:
    if isinstance(shape_names, str):
        shape_names = [shape_names]
    shapes = shapes.query(f"{column} in {shape_names}")
shapes = shapes.set_index(column)

cutout_kwargs = snakemake.params.get("cutout_kwargs", {})
cutout = atlite.Cutout(path=snakemake.input.cutout, **cutout_kwargs)

assert cutout.bounds[0] < shapes.bounds.minx.min()
assert cutout.bounds[1] < shapes.bounds.miny.min()
assert cutout.bounds[2] > shapes.bounds.maxx.max()
assert cutout.bounds[3] > shapes.bounds.maxy.max()

wind_kwargs = snakemake.params.get("wind_kwargs", {})
wind = cutout.wind(
    turbine=snakemake.params.wind_turbine,
    add_cutout_windspeed=True,
    per_unit=True,
    shapes=shapes,
    **wind_kwargs,
)
wind.name = "wind-cf-series"
wind_df = wind.to_dataframe().unstack()
wind_df.columns = wind_df.columns.droplevel(0)
wind_df.to_csv(snakemake.output.timeseries)

plot_mean_cf = snakemake.output.get("plot_mean_cf")
if plot_mean_cf:
    mean = wind.mean("time").to_series()
    shapes.plot(column=mean, legend=True, cmap="PuBu")
    plt.savefig(plot_mean_cf)
