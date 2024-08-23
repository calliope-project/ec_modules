"""Wrapper for atlite availability matrixes."""

import atlite
import geopandas as gpd
import matplotlib.pyplot as plt
from atlite.gis import ExclusionContainer
from pyproj import CRS
from shapely import box

shapes = gpd.read_file(snakemake.input.shapefile)
shapes = shapes.set_index(snakemake.params.shapefile_index_column)
assert shapes.crs.is_geographic
assert shapes.index.is_unique

cutout = atlite.Cutout(path=snakemake.input.cutout)
assert shapes.crs == cutout.crs
assert box(*cutout.bounds).covers(box(*shapes.total_bounds))

exclusion_container_kwargs = snakemake.input.get("exclusion_container_kwargs", {})
excluder = ExclusionContainer(**exclusion_container_kwargs)
assert CRS.from_user_input(excluder.crs).is_projected

rasters = snakemake.input.rasters
raster_codes = snakemake.params.raster_codes
raster_kwargs = snakemake.params.raster_kwargs
assert len(raster_codes) == len(raster_codes) == len(raster_kwargs)
assert all([isinstance(i, list) for i in [rasters, raster_codes, raster_kwargs]])
for i, raster in enumerate(rasters):
    excluder.add_raster(raster, codes=raster_codes[i], **raster_kwargs[i])

availability = cutout.availabilitymatrix(
    shapes, excluder, **snakemake.params.get("availability_kwargs", {})
)
availability.name = snakemake.params.availability_name
availability.attrs["units"] = "%"
availability.to_netcdf(snakemake.output.availability_matrix)

plot_raster = snakemake.output.get("plot_availability_raster")
if plot_raster:
    geometries = shapes.geometry.to_crs(excluder.crs)
    fig, ax = plt.subplots()
    excluder.plot_shape_availability(geometries, ax=ax)
    cutout.grid.to_crs(excluder.crs).plot(edgecolor="grey", color="None", ax=ax, ls=":")
    ax.axis("off")
    plt.savefig(plot_raster)

plot_availability = snakemake.output.get("plot_availability_matrix")
if plot_availability:
    fig, ax = plt.subplots()
    availability.sum(shapes.index.name, keep_attrs=True).plot(cmap="Greens")
    shapes.plot(ax=ax, edgecolor="k", color="None")
    cutout.grid.plot(ax=ax, color="None", edgecolor="grey", ls=":")
    plt.savefig(plot_availability)
