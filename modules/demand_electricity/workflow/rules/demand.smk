"""Rules to generate electricity demand time series."""

rule electricity_load_national:
    message: "Preprocess raw electricity load data and retrieve load time series per country."
    input:
        load = rules.download_raw_load.output[0]
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        data_quality_config = config["quality-control"]["load"],
        countries = config["scope"]["spatial"]["countries"]
    output: "results/electricity-demand-national.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/demand/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location on {wildcards.resolution} resolution."
    input:
        units = rules.download_units.output[0],
        demand_per_unit = rules.unzip_potentials.output.demand,
        national_load = rules.electricity_load_national.output[0]
    params:
        scaling_factor = config["scaling_factors"]["power"]
    output: "results/demand-electricity-{resolution}.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/demand/load.py"
