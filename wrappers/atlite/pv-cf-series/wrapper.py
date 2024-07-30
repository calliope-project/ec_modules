"""atlite PV capacity factor wrapper."""

import atlite
import geopandas as gpd
from matplotlib import pyplot as plt

CUTOUT_OPTIONAL_PARAMS = {"cutout_chunks"}
PV_OPTIONAL_PARAMS = {"pv_tracking", "pv_clearsky_model"}

shapes = gpd.read_file(snakemake.input.shapefile)

column = snakemake.params.shapefile_name_column
shape_names = snakemake.params.get("shapefile_specific_names")
if shape_names:
    if isinstance(shape_names, str):
        shape_names = [shape_names]
    shapes = shapes.query(f"{column} in {shape_names}")
shapes = shapes.set_index(column)

cutout_optional_params = {
    param: snakemake.params.get(param)
    for param in CUTOUT_OPTIONAL_PARAMS
    if param in snakemake.params.keys()
}
cutout = atlite.Cutout(path=snakemake.input.cutout, **cutout_optional_params)

assert cutout.bounds[0] < shapes.bounds.minx.min()
assert cutout.bounds[1] < shapes.bounds.miny.min()
assert cutout.bounds[2] > shapes.bounds.maxx.max()
assert cutout.bounds[3] > shapes.bounds.maxy.max()

pv_optional_params = {
    param: snakemake.params.get(param)
    for param in PV_OPTIONAL_PARAMS
    if param in snakemake.params.keys()
}
pv = cutout.pv(
    panel=snakemake.params.pv_panel,
    orientation=snakemake.params.pv_orientation,
    per_unit=True,
    shapes=shapes,
    **pv_optional_params,
)
pv.name = "pv-cf-series"
pv_df = pv.to_dataframe().unstack()
pv_df.columns = pv_df.columns.droplevel(0)
pv_df.to_csv(snakemake.output.timeseries)

plot_mean_cf = snakemake.output.get("plot_mean_cf")
if plot_mean_cf:
    mean = pv.mean("time").to_series()
    shapes.plot(column=mean, legend=True, cmap="inferno")
    plt.savefig(plot_mean_cf)
