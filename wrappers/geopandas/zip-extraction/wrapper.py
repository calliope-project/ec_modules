"""Extract shapes from a zip file."""

import geopandas as gpd
import matplotlib.pyplot as plt

shapes = gpd.read_file(f"{snakemake.input.zipfile}!{snakemake.params.zip_filepath}")
shapes.to_file(snakemake.output.shapefile, **snakemake.params.get("to_file_kwargs", {}))

plot_shapefile = snakemake.output.get("plot_shapefile")
if plot_shapefile:
    shapes.plot()
    plt.savefig(plot_shapefile)
