"""Rules to generate hydro electricity capacities and time series."""

import yaml


# TODO: should improve the code to automatically determine countries
# and bounds based on the supplied GeoJSON file, and all intermediate
# build artefacts resolution-dependent or unit_id-dependent rather
# than bounds-dependent
internal_config = yaml.safe_load("""
data-sources:
    hydro-generation: https://zenodo.org/record/5797549/files/hydro-generation.csv?download=1
    national-phs-storage-capacities: https://zenodo.org/record/5797549/files/pumped-hydro-storage-capacities-gwh.csv?download=1
    hydro-stations: https://zenodo.org/record/5215920/files/energy-modelling-toolkit/hydro-power-database-v10.zip?download=1

capacity-factors:
    min: 0.001
    max: 10 # 0.001 -> 10 leads to a numerical range of 1e5 (hourly resolution)

quality-control:
    hydro:
        scale-phs-according-to-geth-et-al: false
        station-nearest-basin-max-km: 1
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
    output: "results/downloads/hydro-generation.csv"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule download_pumped_hydro_data:
    message: "Download database of pumped hydro storage capacity data."
    params: url = internal_config["data-sources"]["national-phs-storage-capacities"]
    output: "results/downloads/pumped-hydro-storage-capacities-gwh.csv"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule powerplants_download:
    message: "Download the JRC hydro-power-database."
    params: url = internal_config["data-sources"]["hydro-stations"]
    output: temp("results/downloads/powerplants.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule powerplants_unzip:
    message: "Unzip the JRC hydro-power-database."
    input: rules.powerplants_download.output
    output: "results/downloads/jrc-hydro-power-plant-database.csv"
    shadow: "full"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        """
        unzip -j {input} "**/jrc-hydro-power-plant-database.csv" -d results/downloads
        """


rule preprocess_powerplants:
    message: "Preprocess hydro stations."
    input:
        stations = rules.powerplants_unzip.output[0],
        basins = "results/basins/preprocessed_shape_eu.gpkg",
        phs_storage_capacities = rules.download_pumped_hydro_data.output[0]
    params:
        buffer_size_m = internal_config["quality-control"]["hydro"]["station-nearest-basin-max-km"] * 1000,
        countries = internal_config["scope"]["spatial"]["countries"],
        scale_phs = internal_config["quality-control"]["hydro"]["scale-phs-according-to-geth-et-al"]
    output: "results/jrc-hydro-power-plant-database-preprocessed.csv"
    conda: "../envs/hydro.yaml"
    script: "../scripts/preprocess_powerplants.py"


rule hydro_capacities:
    message: "Determine hydro capacities on {wildcards.resolution} resolution."
    input:
        locations = rules.shapes.output[0],
        plants = rules.preprocess_powerplants.output[0]
    output:
        supply = "results/{resolution}/supply-hydro.csv",
        storage = "results/{resolution}/storage-hydro.csv"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro_capacities.py"


rule prepare_runoff_cutout:
    message: "Prepare atlite cutout with runoff data."
    input:
        shapefile = "results/{resolution}/shapes.geojson"
    output:
        cutout = "results/{resolution}/{first_year}-{final_year}/cutout.nc",
        plot_cutout = "results/{resolution}/{first_year}-{final_year}/plots/plot_cutout.png"
    params:
        time = lambda wc: slice(f"{int(wc.first_year)- 1}-01", f"{wc.final_year}-12"),
        features = ["runoff"],
        offset_degrees = 1,
        module = ["era5"],
    wrapper: "v0.0.2/wrappers/atlite/cutout-prepare"


rule inflow_m3:
    message: "Determine water inflow time series for all hydro electricity between the years {wildcards.first_year} and {wildcards.final_year}."
    input:
        stations = rules.preprocess_powerplants.output[0],
        basins = "results/basins/preprocessed_shape_eu.gpkg",
        runoff = "results/{resolution}/{first_year}-{final_year}/cutout.nc"
    output: "results/{resolution}/{first_year}-{final_year}/hydro-electricity-water-inflow.nc"
    conda: "../envs/hydro.yaml"
    resources:
        runtime = 100
    script: "../scripts/inflow_m3.py"


rule inflow_mwh:
    message: "Determine energy inflow time series for all hydro electricity between the years {wildcards.first_year} and {wildcards.final_year}."
    input:
        stations = "results/{resolution}/{first_year}-{final_year}/hydro-electricity-water-inflow.nc",
        generation = rules.download_hydro_generation_data.output[0]
    params:
        max_capacity_factor = internal_config["capacity-factors"]["max"]
    output: "results/{resolution}/{first_year}-{final_year}/hydro-electricity-with-energy-inflow.nc"
    conda: "../envs/hydro.yaml"
    resources:
        runtime = 100,
        memory = 40000
    script: "../scripts/inflow_mwh.py"


rule capacity_factors_hydro:
    message: "Generate capacityfactor time series for hydro electricity on {wildcards.resolution} resolution."
    input:
        capacities = rules.hydro_capacities.output.supply,
        stations = "results/{resolution}/{first_year}-{final_year}/hydro-electricity-with-energy-inflow.nc",
        locations = rules.shapes.output[0]
    params:
        threshold = internal_config["capacity-factors"]["min"]
    output:
        ror = "results/{resolution}/{first_year}-{final_year}/capacity-factors-run-of-river.csv",
        reservoir = "results/{resolution}/{first_year}-{final_year}/capacity-factors-reservoir.csv"
    resources:
        runtime = 100
    conda: "../envs/hydro.yaml"
    script: "../scripts/capacityfactors_hydro.py"
