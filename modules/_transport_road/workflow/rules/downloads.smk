"""Rules to download necessary files before running the workflow."""

rule download_units:
    message: "Download spatial zones."
    params:
        url = config["data-sources"]["spatial-zones"],
    output: "results/downloads/units.geojson"
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"


rule download_transport_timeseries:
    message: "Get EV data from RAMP"
    params:
        url = internal["data-sources"]["controlled-ev-profiles"]
    conda: "../envs/shell.yaml"
    output: temp("results/downloads/ramp-ev-{dataset}.csv.gz")
    wildcard_constraints:
        dataset = "consumption-profiles|plugin-profiles"
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


# TODO: move into rule download_transport_timeseries once PR 356 is merged
rule download_uncontrolled_timeseries:
    message: "Get EV uncontrolled charging data from RAMP"
    params:
        url = internal["data-sources"]["uncontrolled-ev-profiles"]
    conda: "../envs/shell.yaml"
    output: temp("results/downloads/ramp-ev-uncontrolled-charging-profiles.csv.gz")
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


rule download_jrc_transport_road:
    message: "Get tidy JRC data for road transport."
    params:
        distance_url = internal["data-sources"]["jrc-idees"]["road-distance"],
        energy_url = internal["data-sources"]["jrc-idees"]["road-energy"],
        vehicles_url = internal["data-sources"]["jrc-idees"]["road-vehicles"]
    conda: "../envs/shell.yaml"
    output:
        distance = "results/downloads/jrc-idees-transport-road-distance.csv",
        energy = "results/downloads/jrc-idees-transport-road-energy.csv",
        vehicles = "results/downloads/jrc-idees-transport-road-vehicles.csv"
    localrule: True
    shell:
        """
        curl -sSLo '{output.distance}' '{params.distance_url}'
        curl -sSLo '{output.energy}' '{params.energy_url}'
        curl -sSLo '{output.vehicles}' '{params.vehicles_url}'
        """


rule download_annual_energy_balances:
    message: "Get tidy annual energy balance data for European countries."
    params:
        url = internal["data-sources"]["annual-energy-balances"]
    conda: "../envs/shell.yaml"
    output:
        balances = "results/downloads/annual-energy-balances.csv"
    localrule: True
    shell: "curl -sSLo '{output.balances}' '{params.url}'"
