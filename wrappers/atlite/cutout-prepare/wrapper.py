"""Wrapper for atlite cutout preparation."""

import atlite
import geopandas as gpd
from matplotlib import pyplot as plt

OPTIONAL_INPUTS = ["gebco_path"]

cutout_kwargs = snakemake.params.get("cutout_kwargs", {})
for optional_input in OPTIONAL_INPUTS:
    value = snakemake.input.get(optional_input)
    if value:
        cutout_kwargs[optional_input] = value

shapes = gpd.read_file(snakemake.input.shapefile)
assert shapes.crs.is_geographic

offset = snakemake.params.get("offset_degrees", 0)
assert offset >= 0

cutout = atlite.Cutout(
    path=snakemake.output.cutout,
    module=snakemake.params.module,
    x=slice(shapes.total_bounds[0] - offset, shapes.total_bounds[2] + offset),
    y=slice(shapes.total_bounds[1] - offset, shapes.total_bounds[3] + offset),
    time=snakemake.params.time,
    **cutout_kwargs,
)
cutout.prepare(
    features=snakemake.params.features, **snakemake.params.get("prepare_kwargs", {})
)

plot_cutout = snakemake.output.get("plot_cutout")
if plot_cutout:
    plt.rc("figure", figsize=[10, 7])
    fig, ax = plt.subplots()
    shapes.plot(ax=ax)
    cutout.grid.plot(ax=ax, edgecolor="grey", color="None")
    plt.savefig(plot_cutout)
