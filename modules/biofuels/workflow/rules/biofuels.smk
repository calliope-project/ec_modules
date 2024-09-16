"""Rules related to biofuels."""

rule preprocess_biofuel_potentials_and_cost:
    message: "Extract national potentials and cost from raw biofuel data."
    input:
        potentials_and_costs = rules.download_biofuel_potentials_and_costs.output[0]
    params:
        feedstocks = {
            feedstock["id"]: name
            for name, feedstock in internal["feedstocks"].items()
            if feedstock["include"]
        }
    output:
        potentials = "results/raw_biofuel_potentials.csv",
        costs = "results/raw_biofuel_costs.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/extract.py"


rule biofuels:
    message: "Determine biofuels potential on {wildcards.resolution} resolution for scenario {wildcards.scenario}."
    input:
        units = "resources/user/shapes_{resolution}.geojson",
        land_cover = rules.unzip_potentials.output.land_cover,
        population = rules.unzip_potentials.output.population,
        national_potentials = rules.preprocess_biofuel_potentials_and_cost.output.potentials,
        costs = rules.preprocess_biofuel_potentials_and_cost.output.costs
    params:
        potential_year = internal["potential-year"],
        cost_year = internal["cost-year"],
        proxies = {
            name: feedstock["proxy"]
            for name, feedstock in internal["feedstocks"].items()
            if feedstock["include"]
        }
    output:
        potentials = "results/{resolution}/{scenario}/potential_mwh_per_year.csv",
        costs = "results/{resolution}/{scenario}/costs_eur_per_mwh.csv" # not actually resolution dependent
    conda: "../envs/default.yaml"
    wildcard_constraints:
        scenario = "low|medium|high",
        resolution = "ehighways|national|regional|continental"
    script: "../scripts/allocate.py"

rule plot:
    input:
        shapes = "resources/user/shapes_{resolution}.geojson",
        potentials = rules.biofuels.output.potentials,
        costs = rules.biofuels.output.costs
    output: "results/{resolution}/{scenario}/potentials.png"
    conda: "../envs/default.yaml"
    script: "../scripts/plot.py"
