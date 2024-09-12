"""Rules to generate electricity demand time series."""

rule electricity_load_national:
    message: "Preprocess raw electricity load data and retrieve load time series per country."
    input:
        load = rules.download_raw_load.output[0]
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        data_quality_config = internal["quality-control"]["load"],
        countries = internal["scope"]["spatial"]["countries"]
    output: "results/electricity-demand-national.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location on {wildcards.resolution} resolution."
    input:
        units = "resources/user/shapes_{resolution}.geojson",
        demand_per_unit = rules.unzip_potentials.output.demand,
        national_load = rules.electricity_load_national.output[0]
    params:
        scaling_factor = internal["scaling_factors"]["power"]
    output: "results/demand-electricity-{resolution}.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/load.py"

rule plot:
    input: rules.electricity_load.output
    output: "results/plot_{resolution}.png"
    conda: "../envs/geo.yaml"
    script: "../scripts/plot.py"
