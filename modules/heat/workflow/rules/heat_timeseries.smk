# The timeseries data comes from 4 different types of data: when2heat, gridded weather data, population, and geounits.

# when2heat
rule download_when2heat_params:
    message: "Get parameters for heat demand profiles from the When2Heat project repository"
    output: directory("results/downloads/when2heat")
    params:
        url = lambda wildcards: internal["data-sources"]["when2heat-params"].format(dataset=
            "{" + ",".join(["daily_demand.csv", "hourly_factors_COM.csv", "hourly_factors_MFH.csv", "hourly_factors_SFH.csv"]) + "}"
        )
    conda: "../envs/shell.yaml"
    shell: "mkdir -p {output} && curl -f -sSLo '{output}/#1' '{params.url}'"


# gridded weather data
rule download_gridded_weather_data:
    message: "Download gridded {wildcards.data_var} data"
    params: dataset_url = internal["data-sources"]["gridded-weather-data"]
    output: "results/downloads/gridded-weather/{data_var}.nc"
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -f -sSLo {output} '{params.dataset_url}/files/{wildcards.data_var}.nc'"


# population data
rule download_raw_population_zipped:
    message: "Download population data."
    output: temp("results/downloads/raw-population-data.zip")
    params: url = internal["data-sources"]["population"]
    conda: "../envs/shell.yaml"
    shell: "curl -sSLo {output} '{params.url}'"


rule raw_population_unzipped:
    message: "Extract population data TIF."
    input: rules.download_raw_population_zipped.output
    shadow: "minimal"
    output: "results/downloads/JRC_1K_POP_2018.tif"
    conda: "../envs/shell.yaml"
    shell: "unzip -p '{input}' 'JRC_1K_POP_2018.tif' > '{output}'"


# synthesizing
rule unscaled_heat_profiles:
    message: "Generate gridded heat demand profile shapes from weather and population data"
    input:
        wind_speed = "results/downloads/gridded-weather/wind10m.nc",
        temperature = "results/downloads/gridded-weather/temperature.nc",
        when2heat = rules.download_when2heat_params.output[0]
    params:
        first_year = config["temporal-scope"]["first-year"],
        final_year = config["temporal-scope"]["final-year"],
    conda: "../envs/default.yaml"
    output: "results/hourly_unscaled_heat_demand.nc"
    script: "../scripts/unscaled_heat_profiles.py"

rule population_per_weather_gridbox:
    message: "Get population information per weather data gridbox"
    input:
        weather_grid = "results/downloads/gridded-weather/grid.nc",
        population = rules.raw_population_unzipped.output[0],
        locations = "results/{shapes}/shapes.geojson"
    params:
        lat_name = "lat",
        lon_name = "lon",
    conda: "../envs/geo.yaml"
    output: "results/{shapes}/population.nc"
    script: "../scripts/population_per_gridbox.py"

rule group_gridded_timeseries_heat_demand:
    message: "Generate heat demand hourly timeseries data from gridded data "
    input:
        gridded_timeseries_data = "results/hourly_unscaled_heat_demand.nc",
        grid_weights = "results/{shapes}/population.nc",
    conda: "../envs/default.yaml"
    threads: 4
    output: temp("results/{shapes}/hourly_unscaled_heat_demand.nc")
    script: "../scripts/group_gridded_timeseries.py"
