"""Timeseries data comes from 4 different types of data:

When2heat, gridded weather data, population, and geounits.
"""

# synthesizing
rule unscaled_heat_profiles:
    message: "Generate gridded heat demand profile shapes from weather and population data"
    input:
        wind_speed = "resources/automatic/gridded-weather/wind10m.nc",
        temperature = "resources/automatic/gridded-weather/temperature.nc",
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
        weather_grid = "resources/automatic/gridded-weather/grid.nc",
        population = rules.unzip_raw_population.output[0],
        locations = "resources/user/{shapes}.geojson"
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


rule heat_demand_final_timeseries:
    message: "Generate heat demand timeseries data from gridded data for '{wildcards.shapes}' resolution."
    input:
        timeseries_data = "results/{shapes}/hourly_unscaled_heat_demand.nc",
        annual_demand = "results/{shapes}/annual-heat-demand-twh.csv",
    conda: "../envs/default.yaml"
    params:
        sfh_mfh_shares = config["parameters"]["heat"]["sfh-mfh-shares"],
        scaling_factor = internal["scaling-factors"]["power"]
    output: "results/{shapes}/timeseries/heat_demand.csv"
    script: "../scripts/heat_demand_final_timeseries.py"
