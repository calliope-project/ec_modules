"""This part only contains heat pump specific rules."""

rule download_heat_pump_characteristics:
    message: "Download manufacturer heat pump data"
    params: url = internal["data-sources"]["heat-pump-characteristics"]
    output: "results/heat-pump-characteristics.nc"
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"


rule heat_pump_cop:
    message: "Generate gridded heat pump coefficient of performance (COP)"
    input:
        temperature_air = "results/downloads/gridded-weather/temperature.nc",
        temperature_ground = "results/downloads/gridded-weather/tsoil5.nc",
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
