"""Download and correct files in the HydroSHEDS-v1 dataset.

https://www.hydrosheds.org/products/hydrobasins
Reference in resources/references/hydroBASINS.bib
"""
# TODO: preset to download global combination, or run to combine all regions.
# - This would  avoid arbitrary "chopped regions", like iceland.
# - Rivers do not care about borders.
# - Would ensure that no regions are missed in case of mismatches between HydroSHEDS and the provided shapes.
HYDROSHEDS_CODES = ["af", "ar", "as", "au", "eu", "gr", "na", "sa"]

rule hydrobasin_download:
    message: "Download HydroSHEDS-HydroBASINS v1.0 zip file for '{wildcards.continent}'."
    conda: "../envs/shell.yaml"
    params:
        prefix = "https://data.hydrosheds.org/file/hydrobasins/standard/hybas_",
        suffix = "_lev01-12_v1c.zip"
    output: temp("results/basins/hydrobasin_{continent}.zip")
    shell:
        "curl -sSLo {output} '{params.prefix}{wildcards.continent}{params.suffix}' "


def zip_filepath(wildcards):
    return f"hybas_{wildcards.continent}_lev{config["basins"]["HydroBASINS-level"]}_v1c.shp"
rule hydrobasin_zip_extract:
    message: "Extract requested basin resolution level."
    input:
        zipfile = "results/basins/hydrobasin_{continent}.zip"
    output:
        shapefile = "results/basins/raw_shape_{continent}.geojson",
        plot_shapefile = "results/basins/plots/basin_{continent}.png"
    params:
        zip_filepath = zip_filepath,
    wrapper: "v0.0.1/wrappers/geopandas/zip-extraction"

# FIXME: unclear if necessary. Might be distorting results?
rule preprocess_basins:
    message: "Preprocess basins."
    input:
        basins = "results/basins/raw_shape_eu.geojson"
    params:
        x_min = internal_config["scope"]["spatial"]["bounds"]["x_min"],
        x_max = internal_config["scope"]["spatial"]["bounds"]["x_max"],
        y_min = internal_config["scope"]["spatial"]["bounds"]["y_min"],
        y_max = internal_config["scope"]["spatial"]["bounds"]["y_max"]
    output: "results/basins/preprocessed_shape_eu.gpkg"
    conda: "../envs/hydro.yaml"
    script: "../scripts/preprocess_basins.py"
