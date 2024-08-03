"""Wrapper for atlite availability matrixes."""

import atlite
import geopandas as gpd
import matplotlib.pyplot as plt
from atlite.gis import ExclusionContainer
from pyproj import CRS
from shapely import box

shapes = gpd.read_file(snakemake.input.shapefile)
shapes = shapes.set_index(snakemake.params.shapefile_name_column)
assert shapes.crs.is_geographic

cutout = atlite.Cutout(path=snakemake.input.cutout)
assert shapes.crs == cutout.crs
assert box(*cutout.bounds).covers(box(*shapes.total_bounds))

excluder = ExclusionContainer(
    crs=snakemake.params.exclusion_crs, res=snakemake.params.exclusion_resolution
)
assert CRS.from_user_input(excluder.crs).is_projected

rasters = snakemake.input.rasters
raster_codes = snakemake.params.raster_codes
raster_kwargs = snakemake.params.raster_kwargs
assert len(raster_codes) == len(raster_codes) == len(raster_kwargs)
assert all([isinstance(i, list) for i in [rasters, raster_codes, raster_kwargs]])
for i, raster in enumerate(rasters):
    excluder.add_raster(raster, codes=raster_codes[i], **raster_kwargs[i])

availability = cutout.availabilitymatrix(
    shapes, excluder, **snakemake.params.get("availability_matrix_kwargs", {})
)
availability.to_netcdf(snakemake.output.availability_matrix)

plot_shape = snakemake.output.get("plot_availability_shape")
if plot_shape:
    geometries = shapes.geometry.to_crs(excluder.crs)
    fig, ax = plt.subplots()
    excluder.plot_shape_availability(geometries, ax=ax)
    cutout.grid.to_crs(excluder.crs).plot(edgecolor="grey", color="None", ax=ax, ls=":")
    ax.axis("off")
    plt.savefig(plot_shape)

plot_matrix = snakemake.output.get("plot_availability_matrix")
if plot_matrix:
    fig, ax = plt.subplots()
    availability.sum(shapes.index.name).plot(cmap="Greens")
    shapes.plot(ax=ax, edgecolor="k", color="None")
    cutout.grid.plot(ax=ax, color="None", edgecolor="grey", ls=":")
    plt.savefig(plot_matrix)
