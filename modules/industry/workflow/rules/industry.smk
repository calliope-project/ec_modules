"""Processing individual JRC-IDEES industry categories."""

# if "Iron and steel" in config['specific-categories']:
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
        path_output = "results/categories/annual_energy_demand_iron_and_steel.nc"
    script: "../scripts/categories/iron_and_steel.py"


# if "Chemicals Industry" in config['specific-categories']:
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
        energy = "results/categories/annual_energy_demand_chemicals_industry.nc",
        mass = "results/categories/annual_mass_demand_chemicals_industry.nc"
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
    output: "results/categories/annual_energy_demand_combined_categories.nc"
    script: "../scripts/categories/combined_categories.py"


rule business_as_usual_electricity:
    message: "Obtain current electricity consumption."
    conda: "../envs/industry.yaml"
    input:
        energy_balances = "results/annual_energy_balances.csv",
        cat_names = workflow.source_path("../internal/energy_balance_category_names.csv"),
    output: "results/annual_bau_electricity_demand.nc"
    script: "../scripts/bau_electricity.py"


SUFFIXES = [i.lower().replace(" ", "_") for i in config["specific-categories"]]
rule combine_and_scale:
    message: "Combine and scale industry sectors."
    conda: "../envs/industry.yaml"
    input:
        expand("results/categories/annual_energy_demand_{sample}.nc", sample=SUFFIXES),
        rules.combined_categories.output,
        rules.business_as_usual_electricity.output,
    output: "results/annual_demand_aggregated.nc"
    shell: "touch {output}"

# FIXME: SCEC evaluates the resulting industry energy demand using IQR. Add it.
# rule evaluate_industry:

# FIXME: SCEC scales the resulting demand with user inputs. Add it.
# rule scaling:
