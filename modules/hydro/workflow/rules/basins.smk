"""Download, correct, and combine files in the HydroSHEDS-v1 dataset.

https://www.hydrosheds.org/products/hydrobasins
"""

HYDROSHEDS_CODES = ["af", "ar", "as", "au", "eu", "gr", "na", "sa"]

# Should be in config
HYDROSHEDS_LEVELS = ['01','02','03','04','05','06','07','08','09','10','11','12']

rule hydrobasin_download:
    message: "Download HydroSHEDS-HydroBASINS v1.0 zip file for '{wildcards.region}'."
    conda: "../envs/shell.yaml"
    params:
        prefix = "https://data.hydrosheds.org/file/hydrobasins/standard/hybas_",
        suffix = "_lev01-12_v1c.zip"
    output: "results/basins/hydrobasin_{region}.zip"
    shell:
        "curl -sSLo {output} '{params.prefix}{wildcards.region}{params.suffix}' "

# # No need to unzip, let geopandas take care of it!
# rule hydrobasin_extract:
#     message: "Unzip basins database."
#     input: "results/basins/hydrobasin_{region}.zip"
#     output: "results/basins/{region}/hybas_{region}_lev07_v1c.shp"
#     conda: "../envs/shell.yaml"
#     localrule: True
#     shell: "unzip {input} -d ./results/basins/'{wildcards.region}'/"
