# We recommend adding rules that download necessary files here.

if config["use_default_resources"]:
    rule download_units:
        message:
            "Download spatial zones."
        params:
            url=internal["default_resources"]["spatial_units"],
        output:
            "resources/customizable_resources/units_{resolution}.geojson",
        conda:
            "../envs/shell.yaml"
        localrule: True
        shell:
            "curl -sSLo {output} '{params.url}'"


rule download_eez:
    message:
        "Download Exclusive Economic Zones as zip"
    output:
        protected("resources/eez.gpkg.zip"),
    params:
        url=internal["resources"]["eez"],
    conda:
        "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule download_capacity_factors_wind_and_solar:
    message:
        "Download resources/capacityfactors/{wildcards.filename}."
    params:
        url=lambda wildcards: internal["resources"]["capacity-factors"].format(
            filename=wildcards.filename
        ),
    output:
        protected("resources/capacityfactors/{filename}"),
    conda:
        "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"
