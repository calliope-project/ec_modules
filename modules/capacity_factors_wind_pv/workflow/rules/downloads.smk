# We recommend adding rules that download necessary files here.

if config["use_default_user_resources"]:
    rule user_input_shapes:
        message:
            "Download spatial zones."
        params:
            url=internal["default_user_resources"]["spatial_units"],
        output:
            "resources/user/shapes_ehighways.geojson",
        conda:
            "../envs/shell.yaml"
        localrule: True
        shell:
            "curl -sSLo {output} '{params.url}'"


rule download_eez:
    message:
        "Download Exclusive Economic Zones as zip"
    output:
        protected("resources/automatic/eez.gpkg.zip"),
    params:
        url=internal["resources"]["eez"],
    conda:
        "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule download_capacity_factors_wind_and_solar:
    message:
        "Download capacity factors: {wildcards.filename}."
    params:
        url=lambda wildcards: internal["resources"]["capacity-factors"].format(
            filename=wildcards.filename
        ),
    output:
        protected("resources/automatic/capacityfactors/{filename}"),
    conda:
        "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"
