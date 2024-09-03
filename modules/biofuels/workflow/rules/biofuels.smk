"""Rules related to biofuels."""



rule preprocess_biofuel_potentials_and_cost:
    message: "Extract national potentials and cost from raw biofuel data."
    input:
        potentials_and_costs = rules.download_biofuel_potentials_and_costs.output[0]
    params:
        feedstocks = {
            feedstock["id"]: name
            for name, feedstock in internal["parameters"]["jrc-biofuel"]["feedstocks"].items()
            if feedstock["include"]
        }
    output:
        potentials = "results/raw-biofuel-potentials.csv",
        costs = "results/raw-biofuel-costs.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/extract.py"


rule biofuels:
    message: "Determine biofuels potential on {wildcards.resolution} resolution for scenario {wildcards.scenario}."
    input:
        units = "resources/customisable/spatial_units.geojson",
        land_cover = rules.unzip_potentials.output.land_cover,
        population = rules.unzip_potentials.output.population,
        national_potentials = rules.preprocess_biofuel_potentials_and_cost.output.potentials,
        costs = rules.preprocess_biofuel_potentials_and_cost.output.costs
    params:
        potential_year = internal["parameters"]["jrc-biofuel"]["potential-year"],
        cost_year = internal["parameters"]["jrc-biofuel"]["cost-year"],
        proxies = {
            name: feedstock["proxy"]
            for name, feedstock in internal["parameters"]["jrc-biofuel"]["feedstocks"].items()
            if feedstock["include"]
        }
    output:
        potentials = "results/{resolution}/{scenario}/potential-mwh-per-year.csv",
        costs = "results/{resolution}/{scenario}/costs-eur-per-mwh.csv" # not actually resolution dependent
    conda: "../envs/default.yaml"
    wildcard_constraints:
        scenario = "low|medium|high"
    script: "../scripts/allocate.py"
