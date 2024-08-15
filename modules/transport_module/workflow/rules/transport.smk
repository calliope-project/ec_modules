"""Rules to process transport sector data."""

rule download_transport_timeseries:
    message: "Get EV data from RAMP"
    params:
        url = config["module-data-sources"]["controlled-ev-profiles"]
    conda: "../envs/shell.yaml"
    output: temp("results/automatic/ramp-ev-{dataset}.csv.gz")
    wildcard_constraints:
        dataset = "consumption-profiles|plugin-profiles"
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


rule download_uncontrolled_timeseries:
    # TODO: move into rule download_transport_timeseries once PR 356 is merged
    message: "Get EV uncontrolled charging data from RAMP"
    params:
        url = config["module-data-sources"]["uncontrolled-ev-profiles"]
    conda: "../envs/shell.yaml"
    output: temp("results/automatic/ramp-ev-uncontrolled-charging-profiles.csv.gz")
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


rule annual_transport_demand:
    message: "Calculate future transport energy demand based on JRC IDEES"
    input:
        energy_balances = "resources/annual-energy-balances.csv",
        jrc_road_energy = "resources/jrc-idees/processed-road-energy.csv",
        jrc_road_distance = "resources/jrc-idees/processed-road-distance.csv",
    params:
        fill_missing_values = config["data-pre-processing"]["fill-missing-values"]["jrc-idees"],
        efficiency_quantile = config["parameters"]["transport"]["future-vehicle-efficiency-percentile"],
        uncontrolled_charging_share = config["parameters"]["transport"]["uncontrolled-ev-charging-share"],
    conda: "../envs/default.yaml"
    output:
        road_distance_controlled = temp("results/automatic/annual-road-transport-distance-demand-controlled.csv"),
        road_distance_uncontrolled = temp("results/automatic/annual-road-transport-distance-demand-uncontrolled.csv"),
        road_distance_historically_electrified = temp("results/automatic/annual-road-transport-distance-demand-historic-electrification.csv")
    script: "../scripts/annual_transport_demand.py"


rule create_controlled_road_transport_annual_demand_and_installed_capacities:
    message: "Create annual demand for controlled charging and corresponding charging potentials at a given resolution"
    input:
        annual_controlled_demand = "results/automatic/annual-road-transport-distance-demand-controlled.csv",
        ev_vehicle_number = "resources/jrc-idees/processed-road-vehicles.csv",
        jrc_road_distance = "resources/jrc-idees/processed-road-distance.csv",
        locations = "results/automatic/units.csv",
        populations = "resources/population/population_{resolution}.csv".format(
            resolution=config["module-data-sources"]["resolution"],
        ),  #FIXME: temporary hack
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        power_scaling_factor = config["scaling-factors"]["power"],
        transport_scaling_factor = config["scaling-factors"]["transport"],
        battery_sizes = config["parameters"]["transport"]["ev-battery-sizes"],
        conversion_factors = config["parameters"]["transport"]["road-transport-conversion-factors"],
        countries = config["scope"]["spatial"]["countries"],
        fill_missing_values = config["data-pre-processing"]["fill-missing-values"]["jrc-idees"],
        resolution = config["module-data-sources"]["resolution"] #FIXME: temporary hack
    conda: "../envs/default.yaml"
    output:
        main = "results/electrified-transport.csv",
    script: "../scripts/road_transport_controlled_charging.py"


rule create_controlled_ev_charging_parameters:
    message: "Create timeseries parameters {wildcards.dataset_name} for controlled EV charging at a given resolution"
    input:
        ev_profiles = lambda wildcards: "results/automatic/ramp-ev-consumption-profiles.csv.gz" if "demand" in wildcards.dataset_name else f"results/automatic/ramp-ev-{wildcards.dataset_name}.csv.gz",
        locations = "results/automatic/units.csv",
        populations = "resources/population/population_{resolution}.csv".format(
            resolution=config["module-data-sources"]["resolution"],
        ),  #FIXME: temporary hack
    params:
        demand_range = config["parameters"]["transport"]["monthly-demand-bound-fraction"],
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        country_neighbour_dict = config["data-pre-processing"]["fill-missing-values"]["ramp"],
        countries = config["scope"]["spatial"]["countries"],
        resolution = config["module-data-sources"]["resolution"] #FIXME: temporary hack
    wildcard_constraints:
        dataset_name = "demand-shape-equals|demand-shape-max|demand-shape-min|plugin-profiles"
    conda: "../envs/default.yaml"
    output: "results/timeseries/{dataset_name}-ev.csv"
    script: "../scripts/road_transport_controlled_constraints.py"


rule create_uncontrolled_road_transport_timeseries:
    message: "Create timeseries for road transport demand (uncontrolled charging)"
    input:
        annual_data = "results/automatic/annual-road-transport-distance-demand-uncontrolled.csv",
        timeseries = "results/automatic/ramp-ev-uncontrolled-charging-profiles.csv.gz"
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        power_scaling_factor = config["scaling-factors"]["power"],
        conversion_factor = lambda wildcards: config["parameters"]["transport"]["road-transport-conversion-factors"][wildcards.vehicle_type],
        historic = False,
        countries = config["scope"]["spatial"]["countries"],
        country_neighbour_dict = config["data-pre-processing"]["fill-missing-values"]["ramp"]
    conda: "../envs/default.yaml"
    wildcard_constraints:
        vehicle_type = "light-duty-vehicles|heavy-duty-vehicles|coaches-and-buses|passenger-cars|motorcycles"
    output:
        main = "results/timeseries/timeseries-uncontrolled-{vehicle_type}.csv",
    script: "../scripts/road_transport_timeseries.py"


use rule create_uncontrolled_road_transport_timeseries as create_uncontrolled_road_transport_timeseries_historic_electrification with:
    message: "Create timeseries for historic electrified road transport demand (uncontrolled charging)"
    input:
        annual_data = "results/automatic/annual-road-transport-distance-demand-historic-electrification.csv",
        timeseries = "results/automatic/ramp-ev-uncontrolled-charging-profiles.csv.gz"
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        power_scaling_factor = config["scaling-factors"]["power"],
        conversion_factor = lambda wildcards: config["parameters"]["transport"]["road-transport-conversion-factors"][wildcards.vehicle_type],
        historic = True,
        countries = config["scope"]["spatial"]["countries"],
        country_neighbour_dict = config["data-pre-processing"]["fill-missing-values"]["ramp"],
    output:
        "results/timeseries/timeseries-uncontrolled-{vehicle_type}-historic-electrification.csv"

rule aggregate_timeseries: # TODO consider merge with other rules, as this is tiny atm
    message: "Aggregates uncontrolled charging timeseries for electrified road transport transport for a given spatial resolution"
    input:
        time_series = [
            f'results/timeseries/timeseries-uncontrolled-{vehicle_type}.csv'
            for vehicle_type in config["parameters"]["transport"]["road-transport-conversion-factors"].keys()
        ],
        locations = "results/automatic/units.csv",
        populations = "resources/population/population_{resolution}.csv".format(
            resolution=config["module-data-sources"]["resolution"],
        ),  #FIXME: temporary hack
    params:
        resolution = config["module-data-sources"]["resolution"], #FIXME: temporary hack
    conda: "../envs/default.yaml"
    output:
        "results/timeseries/uncontrolled-electrified-road-transport.csv",
    script: "../scripts/aggregate_timeseries.py"


use rule aggregate_timeseries as aggregate_timeseries_historic_electrified with:
    message: "Aggregates uncontrolled charging timeseries for historically electrified road transport for a given spatial resolution"
    input:
        time_series = (
            "results/timeseries/timeseries-uncontrolled-light-duty-vehicles-historic-electrification.csv",
            "results/timeseries/timeseries-uncontrolled-coaches-and-buses-historic-electrification.csv",
            "results/timeseries/timeseries-uncontrolled-passenger-cars-historic-electrification.csv",),
        locations = "results/automatic/units.csv",
        populations = "resources/population/population_{resolution}.csv".format(
            resolution=config["module-data-sources"]["resolution"],
        ),  #FIXME: temporary hack
    params:
        resolution = config["module-data-sources"]["resolution"], #FIXME: temporary hack
    output:
        "results/timeseries/uncontrolled-road-transport-historic-electrification.csv"
