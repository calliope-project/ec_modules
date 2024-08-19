"""Rules to generate hydro electricity capacities and time series."""

# TODO: should improve the code to automatically determine countries
# and bounds based on the supplied GeoJSON file, and all intermediate
# build artefacts resolution-dependent or unit_id-dependent rather
# than bounds-dependent

rule download_hydro_generation_data:
    message: "Download database of historical hydro power generation."
    params: url = internal["data-sources"]["hydro-generation"]
    output: "results/downloads/hydro-generation.csv"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule download_pumped_hydro_storage_capacity:
    message: "Download database of pumped hydro storage capacity data."
    params: url = internal["data-sources"]["national-phs-storage-capacities"]
    output:
        database = "results/downloads/pumped-hydro-storage-capacities-gwh.csv"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.url}'"


rule powerplants_download:
    message: "Download the JRC hydro-power-database."
    params:
        prefix = "https://github.com//energy-modelling-toolkit/hydro-power-database/raw/",
        version = config["hydro_power_database_version"],
        suffix = "/data/jrc-hydro-power-plant-database.csv"
    output:
        database = "results/downloads/jrc-hydro-power-plant-database.csv"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        "curl -sSLo {output} '{params.prefix}{params.version}{params.suffix}'"

M_TO_KM = 1000
rule preprocess_powerplants:
    message: "Preprocess hydro stations."
    input:
        stations = rules.powerplants_download.output.database,
        basins = "results/basins/preprocessed_shape_eu.gpkg",
        phs_storage_capacities = rules.download_pumped_hydro_storage_capacity.output.database
    params:
        buffer_size_m = internal["quality-control"]["hydro"]["station-nearest-basin-max-km"] * M_TO_KM,
        countries = internal["scope"]["spatial"]["countries"],
        scale_phs = internal["quality-control"]["hydro"]["scale-phs-according-to-geth-et-al"]
    output: "results/jrc-hydro-power-plant-database-preprocessed.csv"
    conda: "../envs/hydro.yaml"
    script: "../scripts/preprocess_powerplants.py"


rule hydro_capacities:
    message: "Determine hydro capacities on {wildcards.resolution} resolution."
    input:
        locations = "results/{resolution}/shapes.geojson",
        plants = "results/jrc-hydro-power-plant-database-preprocessed.csv"
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
        stations = "results/jrc-hydro-power-plant-database-preprocessed.csv",
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
        generation = "results/downloads/hydro-generation.csv"
    params:
        max_capacity_factor = internal["capacity-factors"]["max"]
    output: "results/{resolution}/{first_year}-{final_year}/hydro-electricity-with-energy-inflow.nc"
    conda: "../envs/hydro.yaml"
    resources:
        runtime = 100,
        memory = 40000
    script: "../scripts/inflow_mwh.py"


rule capacity_factors_hydro:
    message: "Generate capacityfactor time series for hydro electricity on {wildcards.resolution} resolution."
    input:
        capacities = "results/{resolution}/supply-hydro.csv",
        stations = "results/{resolution}/{first_year}-{final_year}/hydro-electricity-with-energy-inflow.nc",
        locations = "results/{resolution}/shapes.geojson"
    params:
        threshold = internal["capacity-factors"]["min"]
    output:
        ror = "results/{resolution}/{first_year}-{final_year}/capacity-factors-run-of-river.csv",
        reservoir = "results/{resolution}/{first_year}-{final_year}/capacity-factors-reservoir.csv"
    resources:
        runtime = 100
    conda: "../envs/hydro.yaml"
    script: "../scripts/capacityfactors_hydro.py"
