"""This part only contains heat pump specific rules."""

rule heat_pump_cop:
    message: "Generate gridded heat pump coefficient of performance (COP)"
    input:
        temperature_air = "resources/automatic/gridded-weather/temperature.nc",
        temperature_ground = "resources/automatic/gridded-weather/tsoil5.nc",
        heat_pump_characteristics = rules.download_heat_pump_characteristics.output[0]
    params:
        sink_temperature = config["parameters"]["heat-pump"]["sink-temperature"],
        space_heat_sink_shares = config["parameters"]["heat-pump"]["space-heat-sink-shares"],
        correction_factor = config["parameters"]["heat-pump"]["correction-factor"],
        heat_pump_shares = config["parameters"]["heat-pump"]["heat-pump-shares"],
        first_year = config["temporal-scope"]["first-year"],
        final_year = config["temporal-scope"]["final-year"],
    conda: "../envs/default.yaml"
    output: "results/heat-pump-cop.nc"
    script: "../scripts/heat_pump_cop.py"


rule group_gridded_timeseries_hp_cop:
    message: "Generate national heat-pump-cop timeseries data from gridded data "
    input:
        gridded_timeseries_data = "results/heat-pump-cop.nc",
        grid_weights =  "results/{shapes}/population.nc",
    conda: "../envs/default.yaml"
    threads: 4
    output: "results/{shapes}/heat-pump-cop.nc"
    script: "../scripts/group_gridded_timeseries.py"


rule process_heat_pump_timeseries:
    message: "Combine hot water and space heating characteristics to generate a weighted average national heat pump cop `heat` carrier timeseries."
    input:
        timeseries_data = "results/{shapes}/heat-pump-cop.nc",
        annual_demand = "results/{shapes}/annual-heat-demand-twh.csv"
    conda: "../envs/default.yaml"
    output: "results/{shapes}/timeseries/heat_pump_cop.csv"
    script: "../scripts/heat_pump_final_timeseries.py"
