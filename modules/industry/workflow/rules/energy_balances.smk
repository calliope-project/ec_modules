rule annual_energy_balances:
    message: "Combine and standardise energy balance data from eurostat and CHE."
    input:
        energy_balance = "resources/eurostat/energy_balances.tsv.gz",
        ch_energy_balance = "resources/che/energy_balances.xlsx",
        ch_industry_energy_balance = "resources/che/industry_energy_balances.xlsx",
        cat_names = workflow.source_path("../internal/energy_balance_category_names.csv"),
        carrier_names =  workflow.source_path("../internal/energy_balance_carrier_names.csv")
    output: "results/annual_energy_balances.csv"
    params:
        first_year = 2000
    conda: "../envs/default.yaml"
    script: "../scripts/annual_energy_balances.py"
