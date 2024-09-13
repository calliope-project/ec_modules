"""Rules to generate hydro electricity capacities and time series."""

KM_TO_M = 1000
rule preprocess_powerplants:
    message: "Preprocess hydro stations."
    input:
        stations = rules.download_JRC_hydropower_plants.output.database,
        basins = "results/basins/preprocessed_eu.gpkg",
        phs_storage_capacities = rules.download_Geth2015_PHS_capacity.output.database
    params:
        buffer_size_m = config["powerplant_processing"]["station_nearest_basin_max_km"] * KM_TO_M,
        countries = internal["scope"]["countries"],
        scale_phs = config["powerplant_processing"]["scale_phs_according_to_Geth_et_al"]
    output: "results/jrc_hydropower_plant_database_preprocessed.csv"
    conda: "../envs/hydro.yaml"
    script: "../scripts/preprocess_powerplants.py"


rule capacities_per_shape:
    message: "Determine hydro capacities on {wildcards.resolution} resolution."
    input:
        locations = "resources/user/{resolution}.geojson",
        plants = "results/jrc_hydropower_plant_database_preprocessed.csv"
    output:
        supply = "results/shapes/{resolution}/supply_capacity.csv",
        storage = "results/shapes/{resolution}/storage_capacity.csv"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro_capacities.py"


rule inflow_m3:
    message: "Determine volumetric water inflow time series for hydropower for {wildcards.year}."
    input:
        stations = "results/jrc_hydropower_plant_database_preprocessed.csv",
        basins = "results/basins/preprocessed_eu.gpkg",
        runoff = "resources/automatic/shapes/{resolution}/{year}/cutout.nc"
    output: "results/shapes/{resolution}/{year}/hydropower_inflow_m3.nc"
    conda: "../envs/hydro.yaml"
    resources:
        runtime = 100
    script: "../scripts/inflow_m3.py"


rule inflow_mwh:
    message: "Determine energy inflow time series for hydropower for {wildcards.year}."
    input:
        stations = "results/shapes/{resolution}/{year}/hydropower_inflow_m3.nc",
        generation = "resources/automatic/IRENA_hydro_generation.csv"
    params:
        max_capacity_factor = internal["capacity_factors"]["max"]
    output: "results/shapes/{resolution}/{year}/hydropower_inflow_energy.nc"
    conda: "../envs/hydro.yaml"
    resources:
        runtime = 100,
        memory = 40000
    script: "../scripts/inflow_mwh.py"


rule capacity_factors:
    message: "Generate capacityfactor time series for hydro electricity on {wildcards.resolution} resolution."
    input:
        capacities = "results/shapes/{resolution}/supply_capacity.csv",
        stations = "results/shapes/{resolution}/{year}/hydropower_inflow_energy.nc",
        locations = "resources/user/{resolution}.geojson"
    params:
        threshold = internal["capacity_factors"]["min"]
    output:
        ror = "results/shapes/{resolution}/{year}/capacity_factors_ror.csv",
        reservoir = "results/shapes/{resolution}/{year}/capacity_factors_reservoir.csv"
    resources:
        runtime = 100
    conda: "../envs/hydro.yaml"
    script: "../scripts/capacityfactors_hydro.py"
