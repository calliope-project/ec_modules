# We recommend adding rules that download necessary files here.
if config["use_default_customisable_resources"]:
    rule download_units:
        message: "Download spatial zones."
        params:
            url = internal["default_customisable_resources"]["spatial_zones"],
        output: "results/downloads/units.geojson"
        conda: "../envs/shell.yaml"
        localrule: True
        shell: "curl -sSLo {output} '{params.url}'"


rule download_raw_load:
    message: "Download raw load."
    params: url = internal["data_sources"]["load"]
    output: protected("results/downloads/raw-load-data.csv")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sLo {output} '{params.url}'"

rule download_potentials:
    message: "Download potential data."
    params: url = internal["data_sources"]["potentials"]
    output: protected("data/automatic/raw-potentials.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"

rule unzip_potentials:
    message: "Unzip potentials."
    input: rules.download_potentials.output[0]
    output:
        demand = "build/data/{resolution}/demand.csv",
    conda: "../envs/shell.yaml"
    shell: "unzip {input} '{wildcards.resolution}/*' -d build/data"