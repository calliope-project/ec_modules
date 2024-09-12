# We recommend adding rules that download necessary files here.
if config["use_default_user_resources"]:
    rule download_units:
        message: "Download spatial zones."
        params:
            url = internal["resources"]["default_user_shapes"],
        output: "resources/user/shapes_ehighways.geojson"
        conda: "../envs/shell.yaml"
        localrule: True
        shell: "curl -sSLo {output} '{params.url}'"


rule download_raw_load:
    message: "Download raw load."
    params: url = internal["resources"]["load"]
    output: protected("resources/raw_load_data.csv")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sLo {output} '{params.url}'"

rule download_potentials:
    message: "Download potential data."
    params: url = internal["resources"]["potentials"]
    output: protected("resources/raw_potentials.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"

rule unzip_potentials:
    message: "Unzip potentials."
    input: rules.download_potentials.output[0]
    output:
        demand = "resources/{resolution}/demand.csv",
    conda: "../envs/shell.yaml"
    shell: "unzip {input} '{wildcards.resolution}/demand.csv' -d resources"
