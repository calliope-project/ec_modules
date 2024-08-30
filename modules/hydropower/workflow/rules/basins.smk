"""Download and preprocess HydroBASINS v1 data.

https://www.hydrosheds.org/products/hydrobasins
Reference in resources/references/hydroBASINS.bib
"""

rule hydrobasin_zip_extract:
    message: "Extract requested basin resolution level."
    input:
        zipfile = "resources/basins/hydrobasin_{continent}.zip"
    params:
        zip_filepath = lambda wc: f"hybas_{wc.continent}_lev{config["HydroBASINS_level"]}_v1c.shp"
    output:
        shapefile = "resources/basins/raw_{continent}.geojson",
        plot_shapefile = "results/plots/basins/basin_{continent}.png"
    localrule: True
    wrapper: "v0.0.4/wrappers/geopandas/zip-extraction"


rule preprocess_basins:
    message: "Preprocess basins."
    input:
        basins = "resources/basins/raw_eu.geojson"
    output: "results/basins/preprocessed_eu.gpkg"
    conda: "../envs/hydro.yaml"
    script: "../scripts/preprocess_basins.py"
