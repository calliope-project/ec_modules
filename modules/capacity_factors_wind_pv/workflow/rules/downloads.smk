# We recommend adding rules that download necessary files here.

rule download_eez:
    message: "Download Exclusive Economic Zones as zip"
    output: protected("data/automatic/eez.gpkg.zip")
    params: url = config["data-sources"]["eez"]
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"

rule download_capacity_factors_wind_and_solar:
    message: "Download data/automatic/capacityfactors/{wildcards.filename}."
    params: url = lambda wildcards: config["data-sources"]["capacity-factors"].format(filename=wildcards.filename)
    output: protected("data/automatic/capacityfactors/{filename}")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"