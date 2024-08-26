"""Processing individual JRC-IDEES industry categories."""

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
        path_output = "results/industry_categories/annual_demand_iron_and_steel.nc"
    script: "../scripts/categories/iron_and_steel.py"


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
    output: "results/industry_categories/annual_demand_combined_categories.nc"
    script: "../scripts/categories/combined_categories.py"


# SUFFIXES = [i.lower().replace(" ", "_") for i in config["params"]["specific-categories"]]
# rule combine_and_scale:
#     message: "Identify the category scripts to run based on the configuration."
#     conda: "../envs/industry.yaml"
#     input:
#         expand("{path}/annual_demand_{sample}.nc", path=[BUILD_PATH], sample=SUFFIXES),
#         rules.combined_categories.output
#     output: "{BUILD_PATH}/annual_demand_aggregated.nc"
#     shell: "touch {output}"
