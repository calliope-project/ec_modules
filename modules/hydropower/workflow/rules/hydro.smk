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
        locations = "resources/user_input/{resolution}.geojson",
        plants = "results/jrc_hydropower_plant_database_preprocessed.csv"
    output:
        supply = "results/shapes/{resolution}/hydropower_supply_capacity.csv",
        storage = "results/shapes/{resolution}/hydropower_storage_capacity.csv"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro_capacities.py"


rule prepare_ERA5_runoff_cutout:
    message: "Prepare atlite cutout with runoff data for {params.time}."
    input:
        shapefile = "resources/user_input/{resolution}.geojson"
    output:
        cutout = "results/shapes/{resolution}/{year}/cutout.nc",
        plot_cutout = "results/plots/{resolution}/{year}/plot_cutout.png"
    params:
        time = lambda wc: slice(f"{int(wc.year)- config["year_shift"]}-01", f"{wc.year}-12"),
        features = ["runoff"],
        offset_degrees = 0,
        module = ["era5"],
        prepare_kwargs = {"monthly_requests": True, "concurrent_requests": True, "compression": None}
    wrapper: "v0.0.4/wrappers/atlite/cutout"


rule inflow_m3:
    message: "Determine volumetric water inflow time series for hydropower for {wildcards.year}."
    input:
        stations = "results/jrc_hydropower_plant_database_preprocessed.csv",
        basins = "results/basins/preprocessed_eu.gpkg",
        runoff = "results/shapes/{resolution}/{year}/cutout.nc"
    output: "results/shapes/{resolution}/{year}/hydropower_inflow_m3.nc"
    conda: "../envs/hydro.yaml"
    resources:
        runtime = 100
    script: "../scripts/inflow_m3.py"


rule inflow_mwh:
    message: "Determine energy inflow time series for hydropower for {wildcards.year}."
    input:
        stations = "results/shapes/{resolution}/{year}/hydropower_inflow_m3.nc",
        generation = "resources/IRENA_hydro_generation.csv"
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
        capacities = "results/shapes/{resolution}/hydropower_supply_capacity.csv",
        stations = "results/shapes/{resolution}/{year}/hydropower_inflow_energy.nc",
        locations = "resources/user_input/{resolution}.geojson"
    params:
        threshold = internal["capacity_factors"]["min"]
    output:
        ror = "results/shapes/{resolution}/{year}/capacity_factors_ror.csv",
        reservoir = "results/shapes/{resolution}/{year}/capacity_factors_reservoir.csv"
    resources:
        runtime = 100
    conda: "../envs/hydro.yaml"
    script: "../scripts/capacityfactors_hydro.py"
