"""Processing individual JRC-IDEES industry categories."""

if "Iron and steel" in config['specific-categories']:
    rule iron_and_steel:
        message: "Calculate energy demand for the 'Iron and steel' sector in JRC-IDEES."
        conda: "../envs/industry.yaml"
        params:
            config = config["config-iron-and-steel"]
        input:
            energy_balances = "results/annual_energy_balances.csv",
            carrier_names = workflow.source_path("../internal/energy_balance_carrier_names.csv"),
            cat_names = workflow.source_path("../internal/energy_balance_category_names.csv"),
            jrc_industry_energy = "results/processed_industry_energy.nc",
            jrc_industry_production = "results/processed_industry_production.nc",
        output:
            path_output = "results/categories/annual_demand_iron_and_steel.nc"
        script: "../scripts/categories/iron_and_steel.py"


if "Chemicals Industry" in config['specific-categories']:
    rule chemicals_industry:
        message: "Calculate energy demand for the 'Chemicals Industry' sector in JRC-IDEES."
        conda: "../envs/industry.yaml"
        # params:
        #     config = config["params"]["config-chemicals-industry"]
        input:
            energy_balances = "results/annual_energy_balances.csv",
            carrier_names = workflow.source_path("../internal/energy_balance_carrier_names.csv"),
            cat_names = workflow.source_path("../internal/energy_balance_category_names.csv"),
            jrc_industry_energy = "results/processed_industry_energy.nc",
            jrc_industry_production = "results/processed_industry_production.nc",
        output:
            energy = "results/categories/annual_demand_chemicals_industry_energy.nc",
            mass = "results/categories/annual_demand_chemicals_industry_mass.nc"
        script: "../scripts/categories/chemicals_industry.py"


rule combined_categories:
    message: "Calculate energy demand for all other industry sectors in JRC-IDEES."
    conda: "../envs/industry.yaml"
    params:
        specific_categories = config["specific-categories"],
        config = config["config-combined-categories"],
    input:
        energy_balances = "results/annual_energy_balances.csv",
        carrier_names = workflow.source_path("../internal/energy_balance_carrier_names.csv"),
        cat_names = workflow.source_path("../internal/energy_balance_category_names.csv"),
        jrc_industry_energy = "results/processed_industry_energy.nc",
        jrc_industry_production = "results/processed_industry_production.nc",
    output: "results/categories/annual_demand_combined_categories.nc"
    script: "../scripts/categories/combined_categories.py"


rule business_as_usual_electricity:
    message: "Obtain current electricity consumption."
    conda: "../envs/industry.yaml"
    input:
        energy_balances = "results/annual_energy_balances.csv",
        cat_names = workflow.source_path("../internal/energy_balance_category_names.csv"),
    output: "results/annual_bau_electricity_demand.nc"
    script: "../scripts/bau_electricity.py"

