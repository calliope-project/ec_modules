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


rule download_capacity_factors:
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

rule download_potentials:
    message: "Download potentials data."
    params: url = internal["resources"]["potentials"]
    output: protected("resources/automatic/raw-potentials.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"

rule unzip_area_potentials:
    message: "Unzip lang eligibility potentials."
    input: ancient(rules.download_potentials.output[0])
    shadow: "minimal"
    output:
        land_eligibility_km2 = "resources/automatic/{resolution}/{scenario}/areas.csv"
    conda: "../envs/shell.yaml"
    wildcard_constraints:
        scenario = "|".join(internal["scenarios"]),
        resolution = "|".join(internal["resolutions"])
    shell: "unzip -p '{input}' '{wildcards.resolution}/{wildcards.scenario}/areas.csv' > '{output.land_eligibility_km2}'"
