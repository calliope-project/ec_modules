"""Aggregated annual heat demand data comes from 4 sources:

JRC IDEES, eurostat for European countries, Swiss data, and a certain source for population.
Apart from that, there's a reformatting of shapes.geojson file.
"""

# JRC IDEES (subject to change)
JRC_IDEES_SPATIAL_SCOPE = [
    "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL", "ES", "FI", "FR",
    "HR", "HU", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO",
    "SE", "SI", "SK", "UK"
]

rule process_jrc_idees_tertiary:
    message: "Process tertiary heat data from JRC-IDEES"
    input:
        data = expand(
            "resources/automatic/jrc-idees/tertiary_{country_code}.xlsx",
            country_code=JRC_IDEES_SPATIAL_SCOPE
        )
    output: "results/jrc-idees/tertiary_processed.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/jrc-idees-heat.py"


# synthesize eurostat and Swiss data
rule process_annual_energy_balances:
    message: "Get annual energy balances from Eurostat"
    input:
        energy_balance = "resources/automatic/eurostat/energy-balance.tsv.gz",
        ch_energy_balance = "resources/automatic/CHE/energy-balance.xlsx",
        ch_industry_energy_balance = "resources/automatic/CHE/industry-energy-balance.xlsx",
        cat_names = workflow.source_path("../internal/energy-balance-category-names.csv"),
        carrier_names = workflow.source_path("../internal/energy-balance-carrier-names.csv")
    output: temp("results/annual-energy-balances.csv")
    params:
        first_year = 2000
    conda: "../envs/default.yaml"
    script: "../scripts/annual_energy_balance.py"


# TODO: this is not needed. Just use geopandas.
rule extract_shape_data:
    message: "Dataset of units without geo information."
    input:
        units = "resources/user/{shapes}.geojson"
    output: "results/{shapes}/shape_data.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/nogeo.py"


# synthesizing eurostat, Swiss data and JRC IDEES
rule process_annual_heat_demand:
    message: "Calculate national heat demand for household and commercial sectors"
    input:
        hh_end_use = "resources/automatic/eurostat/hh-end-use.tsv.gz",
        ch_end_use = "resources/automatic/CHE/end-use.xlsx",
        energy_balance = rules.process_annual_energy_balances.output[0],
        commercial_demand = "results/jrc-idees/tertiary_processed.csv",
        carrier_names = workflow.source_path("../internal/energy-balance-carrier-names.csv")
    params:
        heat_tech_params = config["parameters"]["heat"]["tech-efficiencies"],
        countries = internal["scope"]["spatial"]["countries"],
        fill_missing_values = internal["data-pre-processing"]["fill-missing-values"]["jrc-idees"]
    conda: "../envs/default.yaml"
    output:
        total_demand = "results/annual-heat-demand-twh.csv",
        electricity = "results/annual-heat-electricity-demand-twh.csv",
    script: "../scripts/annual_heat_demand.py"


# further synthesize with units and population
rule rescale_annual_heat_demand_to_resolution:
    message: "Scale national heat demand for household and commercial sectors"
    input:
        annual_demand = rules.process_annual_heat_demand.output.total_demand,
        electricity = rules.process_annual_heat_demand.output.electricity,
        locations = "results/{shapes}/shape_data.csv",
        populations = "resources/automatic/{shapes}/population_potentials.csv"
    conda: "../envs/default.yaml"
    output:
        total_demand = "results/{shapes}/annual-heat-demand-twh.csv",
        electricity = "results/{shapes}/annual-historic-electrified-heat-demand-twh.csv",
    script: "../scripts/rescale.py"
