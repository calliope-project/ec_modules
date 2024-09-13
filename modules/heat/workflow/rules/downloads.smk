"""Rules to download necessary resources."""

if config["use_default_user_resources"]:
    # Download files for the default example here.
    rule user_input_shapes:
        message: "Download default user input shapefile."
        params:
            uri = lambda wc: internal["resources"]["default_user_shapes"],
        output: "resources/user/national.geojson"
        conda: "../envs/shell.yaml"
        localrule: True
        shell: "curl -sSLo {output} '{params.uri}'"

# Download other resources below this line.
# when2heat
rule download_when2heat_params:
    message: "Get parameters for heat demand profiles from the When2Heat project repository"
    output: directory("resources/automatic/when2heat")
    params:
        url = lambda wildcards: internal["resources"]["when2heat-params"].format(dataset=
            "{" + ",".join(["daily_demand.csv", "hourly_factors_COM.csv", "hourly_factors_MFH.csv", "hourly_factors_SFH.csv"]) + "}"
        )
    conda: "../envs/shell.yaml"
    shell: "mkdir -p {output} && curl -f -sSLo '{output}/#1' '{params.url}'"


# gridded weather data
rule download_gridded_weather_data:
    message: "Download gridded {wildcards.data_var} data"
    params: dataset_url = internal["resources"]["gridded-weather-data"]
    output: "resources/automatic/gridded-weather/{data_var}.nc"
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -f -sSLo {output} '{params.dataset_url}/files/{wildcards.data_var}.nc'"


rule download_raw_population:
    message: "Download population data."
    output: temp("resources/automatic/raw-population-data.zip")
    params: url = internal["resources"]["population"]
    conda: "../envs/shell.yaml"
    shell: "curl -sSLo {output} '{params.url}'"


rule unzip_raw_population:
    message: "Extract population data TIF."
    input: rules.download_raw_population.output
    shadow: "minimal"
    output: "resources/automatic/JRC_1K_POP_2018.tif"
    conda: "../envs/shell.yaml"
    shell: "unzip -p '{input}' 'JRC_1K_POP_2018.tif' > '{output}'"


rule download_jrc_idees:
    message: "Download JRC IDEES zip file for {wildcards.country_code}"
    params: dataset_url = internal["resources"]["jrc-idees"]
    output: temp("resources/automatic/jrc-idees/{country_code}.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -f -sSLo {output} '{params.dataset_url}/JRC-IDEES-2015_All_xlsx_{wildcards.country_code}.zip'"


rule unzip_jrc_idees:
    message: "Unzip JRC-IDEES tertiary sector data for {wildcards.country_code}"
    input:
        country_data = "resources/automatic/jrc-idees/{country_code}.zip",
    output: "resources/automatic/jrc-idees/tertiary_{country_code}.xlsx"
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "unzip -p {input.country_data} JRC-IDEES-2015_Tertiary_{wildcards.country_code}.xlsx > {output}"


rule download_eurostat_energy_data:
    message: "Download {wildcards.dataset} Eurostat data from euro-calliope datasets"
    params:
        url = lambda wc: internal["resources"]["eurostat"][f"{wc.dataset}"]
    wildcard_constraints:
        dataset = "energy-balance|hh-end-use"
    conda: "../envs/shell.yaml"
    output: "resources/automatic/eurostat/{dataset}.tsv.gz"
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


rule download_CHE_energy_data:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: internal["resources"]["CHE"][f"{wildcards.dataset}"]
    output: "resources/automatic/CHE/{dataset}.xlsx"
    conda: "../envs/shell.yaml"
    wildcard_constraints:
        dataset = "energy-balance|industry-energy-balance|end-use"
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


# Tim's population from 'potentials' workflow. For now, only national resolution
rule download_potentials:
    message: "Download potential data."
    params: url = internal["resources"]["potentials"]
    output: temp("resources/automatic/{shapes}/raw-potentials.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"


# TODO: perhaps this should be a user input instead, since it's tied to the shapes?
rule unzip_potentials:
    message: "Unzip potentials."
    input: rules.download_potentials.output[0]
    shadow: "minimal"
    output: "resources/automatic/{shapes}/population_potentials.csv",
    conda: "../envs/shell.yaml"
    wildcard_constraints:
        shapes = "national"
    shell: "unzip -p {input} {wildcards.shapes}/population.csv > {output}"


rule download_heat_pump_characteristics:
    message: "Download manufacturer heat pump data"
    params: url = internal["resources"]["heat-pump-characteristics"]
    output: "resources/automatic/heat-pump-characteristics.nc"
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"
