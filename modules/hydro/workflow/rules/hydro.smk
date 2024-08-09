"""Rules to generate hydro electricity capacities and time series."""

import yaml

internal_config = yaml.safe_load("""
data-sources:
    hydro-generation: https://zenodo.org/record/5797549/files/hydro-generation.csv?download=1
    national-phs-storage-capacities: https://zenodo.org/record/5797549/files/pumped-hydro-storage-capacities-gwh.csv?download=1
    hydro-basins: https://www.dropbox.com/sh/hmpwobbz9qixxpe/AADeU9iCgMd3ZO1KgrFmfWu6a/HydroBASINS/standard/eu/hybas_eu_lev07_v1c.zip?dl=1
    hydro-stations: https://zenodo.org/record/5215920/files/energy-modelling-toolkit/hydro-power-database-v10.zip?download=1

capacity-factors:
    min: 0.001
    max: 10 # 0.001 -> 10 leads to a numerical range of 1e5 (hourly resolution)

quality-control:
    hydro:
        scale-phs-according-to-geth-et-al: false
        station-nearest-basin-max-km: 1

# TODO: should improve the code to automatically determine countries
# and bounds based on the supplied GeoJSON file, and all intermediate
# build artefacts resolution-dependent or unit_id-dependent rather
# than bounds-dependent
scope:
    spatial:
        countries:
            - "Austria"
            - "Belgium"
            - "Bulgaria"
            - "Croatia"
            - "Cyprus"
            - "Czech Republic"
            - "Denmark"
            - "Estonia"
            - "Finland"
            - "France"
            - "Germany"
            - "Greece"
            - "Hungary"
            - "Ireland"
            - "Italy"
            - "Latvia"
            - "Lithuania"
            - "Luxembourg"
            - "Netherlands"
            - "Poland"
            - "Portugal"
            - "Romania"
            - "Slovakia"
            - "Slovenia"
            - "Spain"
            - "Sweden"
            - "United Kingdom"
            - "Albania"
            - "Bosnia and Herzegovina"
            - "Macedonia, Republic of"
            - "Montenegro"
            - "Norway"
            - "Serbia"
            - "Switzerland"
        bounds:
            x_min: -15.8
            x_max: 37
            y_min: 30
            y_max: 75
""")

rule download_hydro_generation_data:
    message: "Download database of historical hydro power generation."
    params: url = internal_config["data-sources"]["hydro-generation"]
    output:
        protected("build/raw-hydro-generation.csv")
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule download_pumped_hydro_data:
    message: "Download database of pumped hydro storage capacity data."
    params: url = internal_config["data-sources"]["national-phs-storage-capacities"]
    output:
        protected("build/raw-pumped-hydro-storage-capacities-gwh.csv")
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule download_runoff_data:
    message: "Create an atlite cutout of Europe consisting of ERA5 runoff data between the years {wildcards.first_year} and {wildcards.final_year}."
    params:
        x_min = internal_config["scope"]["spatial"]["bounds"]["x_min"],
        x_max = internal_config["scope"]["spatial"]["bounds"]["x_max"],
        y_min = internal_config["scope"]["spatial"]["bounds"]["y_min"],
        y_max = internal_config["scope"]["spatial"]["bounds"]["y_max"]
    output:
        protected("build/europe-cutout-{first_year}-{final_year}.nc")
    conda: "../envs/hydro.yaml"
    resources:
        runtime = 240
    script: "../scripts/runoff.py"


rule download_basins_database:
    message: "Download database of hydro basins."
    params: url = internal_config["data-sources"]["hydro-basins"]
    output:
        protected("build/raw-hydro-basins.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule download_stations_database:
    message: "Download database of hydro electricity stations."
    params: url = internal_config["data-sources"]["hydro-stations"]
    output:
        protected("build/raw-hydro-stations.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule basins_database:
    message: "Unzip basins database."
    input: rules.download_basins_database.output
    output: "build/basins/hybas_eu_lev07_v1c.shp"
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "unzip {input} -d ./build/basins/"


rule stations_database:
    message: "Unzip stations database."
    input: rules.download_stations_database.output
    output: "build/jrc-hydro-power-plant-database.csv"
    shadow: "full"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        """
        unzip -j {input} "**/jrc-hydro-power-plant-database.csv" -d build/
        """


rule preprocess_basins:
    message: "Preprocess basins."
    input:
        basins = rules.basins_database.output[0]
    params:
        x_min = internal_config["scope"]["spatial"]["bounds"]["x_min"],
        x_max = internal_config["scope"]["spatial"]["bounds"]["x_max"],
        y_min = internal_config["scope"]["spatial"]["bounds"]["y_min"],
        y_max = internal_config["scope"]["spatial"]["bounds"]["y_max"]
    output: "build/hybas_eu_lev07_v1c.gpkg"
    conda: "../envs/hydro.yaml"
    script: "../scripts/preprocess_basins.py"


rule preprocess_hydro_stations:
    message: "Preprocess hydro stations."
    input:
        stations = rules.stations_database.output[0],
        basins = rules.preprocess_basins.output[0],
        phs_storage_capacities = rules.download_pumped_hydro_data.output[0]
    params:
        buffer_size_m = internal_config["quality-control"]["hydro"]["station-nearest-basin-max-km"] * 1000,
        countries = internal_config["scope"]["spatial"]["countries"],
        scale_phs = internal_config["quality-control"]["hydro"]["scale-phs-according-to-geth-et-al"]
    output: "build/jrc-hydro-power-plant-database-preprocessed.csv"
    conda: "../envs/hydro.yaml"
    script: "../scripts/preprocess_hydro_stations.py"


rule inflow_m3:
    message: "Determine water inflow time series for all hydro electricity between the years {wildcards.first_year} and {wildcards.final_year}."
    input:
        stations = rules.preprocess_hydro_stations.output[0],
        basins = rules.preprocess_basins.output[0],
        runoff = rules.download_runoff_data.output[0]
    output: "build/hydro-electricity-with-water-inflow-{first_year}-{final_year}.nc"
    conda: "../envs/hydro.yaml"
    resources:
        runtime = 100
    script: "../scripts/inflow_m3.py"


rule inflow_mwh:
    message: "Determine energy inflow time series for all hydro electricity between the years {wildcards.first_year} and {wildcards.final_year}."
    input:
        stations = rules.inflow_m3.output[0],
        generation = rules.download_hydro_generation_data.output[0]
    params:
        max_capacity_factor = internal_config["capacity-factors"]["max"]
    output: "build/hydro-electricity-with-energy-inflow-{first_year}-{final_year}.nc"
    conda: "../envs/hydro.yaml"
    resources:
        runtime = 100,
        memory = 40000
    script: "../scripts/inflow_mwh.py"


rule hydro_capacities:
    message: "Determine hydro capacities on {wildcards.resolution} resolution."
    input:
        locations = rules.units.output[0],
        plants = rules.preprocess_hydro_stations.output[0]
    output:
        supply = "results/{resolution}/supply-hydro.csv",
        storage = "results/{resolution}/storage-hydro.csv"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro_capacities.py"


rule capacity_factors_hydro:
    message: "Generate capacityfactor time series for hydro electricity on {wildcards.resolution} resolution."
    input:
        capacities = rules.hydro_capacities.output.supply,
        stations = "build/hydro-electricity-with-energy-inflow-{first_year}-{final_year}.nc".format(
            first_year = config["first-year"],
            final_year = config["final-year"]
        ),
        locations = rules.units.output[0]
    params:
        threshold = internal_config["capacity-factors"]["min"]
    output:
        ror = "results/{resolution}/capacityfactors-hydro-run-of-river-{first_year}-{final_year}.csv",
        reservoir = "results/{resolution}/capacityfactors-hydro-reservoir-{first_year}-{final_year}.csv"
    resources:
        runtime = 100
    conda: "../envs/hydro.yaml"
    script: "../scripts/capacityfactors_hydro.py"
