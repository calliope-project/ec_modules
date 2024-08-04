"""Get a specific natural earth shapefile."""

import cartopy.io.shapereader as shpreader
import geopandas as gpd
import matplotlib.pyplot as plt

shpfilename = shpreader.natural_earth(
    resolution=snakemake.params.resolution,
    category=snakemake.params.category,
    name=snakemake.params.name,
)
regions = gpd.read_file(shpfilename.absolute())
regions.to_file(
    snakemake.output.shapefile, driver=snakemake.params.get("driver", "GeoJSON")
)

plot_shapefile = snakemake.output.get("plot_shapefile")
if plot_shapefile:
    regions.plot()
    plt.savefig(plot_shapefile)
