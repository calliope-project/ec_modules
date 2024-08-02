"""Wrapper for atlite availability matrixes."""

import atlite
import geopandas as gpd
import matplotlib.pyplot as plt
from atlite.gis import ExclusionContainer

shapes = gpd.read_file(snakemake.input.shapefile)
shapes = shapes.set_index(snakemake.params.shapefile_name_column)
cutout = atlite.Cutout(
    path=snakemake.input.cutout, **snakemake.params.get("cutout_kwargs", {})
)

excluder = ExclusionContainer(**snakemake.params.get("exclusion_kwargs", {}))
rasters = snakemake.input.rasters
raster_codes = snakemake.params.raster_codes
raster_kwargs = snakemake.params.raster_kwargs
assert len(raster_codes) == len(raster_codes) == len(raster_kwargs)

if isinstance(rasters, str):
    rasters = [rasters]
for i, raster in enumerate(rasters):
    excluder.add_raster(raster, codes=raster_codes[i], **raster_kwargs[i])
availability = cutout.availabilitymatrix(shapes, excluder)

availability.to_netcdf(snakemake.output.availability_matrix)

if len(snakemake.output) > 1:
    geometries = shapes.geometry.to_crs(excluder.crs)
    plot_shape_availability = snakemake.output.get("plot_shape_availability")
    if plot_shape_availability:
        fig, ax = plt.subplots()
        excluder.plot_shape_availability(geometries, ax=ax)
        cutout.grid.to_crs(excluder.crs).plot(
            edgecolor="grey", color="None", ax=ax, ls=":"
        )
        ax.axis("off")
        plt.savefig(plot_shape_availability)

    plot_availability_matrix = snakemake.output.get("plot_availability_matrix")
    if plot_availability_matrix:
        fig, ax = plt.subplots()
        availability.sum(shapes.index.name).plot(cmap="Greens")
        shapes.plot(ax=ax, edgecolor="k", color="None")
        cutout.grid.plot(ax=ax, color="None", edgecolor="grey", ls=":")
        plt.savefig(plot_availability_matrix)
