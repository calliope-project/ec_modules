"""Rules to download necessary resources."""

if config["use_default_user_input"]:
    rule user_input_shapes:
        message: "Download default national resolution shapes."
        params:
            url = internal["resources"]["default_user_input"]["shapes"],
        output: "resources/user_input/national.geojson"
        conda: "../envs/shell.yaml"
        shell: "curl -sSLo {output} '{params.url}'"


rule download_HydroBASINS:
    message: "Download HydroSHEDS-HydroBASINS v1.0 zip file for '{wildcards.continent}'."
    conda: "../envs/shell.yaml"
    params:
        url = lambda wc: internal["resources"]["HydroBASINS"].format(continent=wc.continent),
    output: temp("resources/basins/hydrobasin_{continent}.zip")
    wildcard_constraints:
        continent = "|".join(["af", "ar", "as", "au", "eu", "gr", "na", "sa"])
    shell:
        "curl -sSLo {output} '{params.url}' "


rule download_IRENA_hydro_generation:
    message: "Download IRENA database of historical hydro power generation."
    params: url = internal["resources"]["IRENA_hydro_generation"]
    output: "resources/IRENA_hydro_generation.csv"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule download_Geth2015_PHS_capacity:
    message: "Download database of pumped hydro storage capacity data."
    params: url = internal["resources"]["Geth2015_PHS_capacity"]
    output:
        database = "resources/Geth_2015_PHS_capacity.csv"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule download_JRC_hydropower_plants:
    message: "Download the JRC hydro-power-database."
    params:
        url = internal["resources"]["JRC_hydropower_plants"].format(
            version=config["JRC_hydropower_plants_version"]
        )
    output:
        database = "resources/jrc_hydropower_plants.csv"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"
