# The aggregated annual heat demand data comes from 4 sources: JRC IDEES, eurostat for European countries, Swiss data, and a certain source for population (don't ask me why).
# Apart from that, there's a reformatting of shapes.geojson file.

# JRC IDEES (subject to change)
JRC_IDEES_SPATIAL_SCOPE = [
    "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL", "ES", "FI", "FR",
    "HR", "HU", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO",
    "SE", "SI", "SK", "UK"
]

rule download_jrc_idees_zipped:
    message: "Download JRC IDEES zip file for {wildcards.country_code}"
    params: url = internal["data-sources"]["jrc-idees"]
    output: "results/jrc-idees/downloads/{country_code}.zip"
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"

rule jrc_idees_unzipped:
    message: "Unzip JRC-IDEES tertiary sector data for {wildcards.country_code}"
    input:
        country_data = "results/jrc-idees/downloads/{country_code}.zip",
    output: temp("results/jrc-idees/tertiary/unprocessed_{country_code}.xlsx")
    conda: "../envs/shell.yaml"
    shell: "unzip -p {input.country_data} JRC-IDEES-2015_Tertiary_{wildcards.country_code}.xlsx > {output}"

rule jrc_idees_tertiary_processed:
    message: "Process tertiary heat data from JRC-IDEES"
    input:
        data = expand(
            "results/jrc-idees/tertiary/unprocessed_{country_code}.xlsx",
            country_code=JRC_IDEES_SPATIAL_SCOPE
        )
    output: "results/jrc-idees/tertiary/processed.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/jrc-idees-heat.py"


# eurostat
rule download_eurostat_energy_data:
    message: "Download {wildcards.dataset} Eurostat data from euro-calliope datasets"
    params:
        url = lambda wildcards: internal["data-sources"][f"eurostat-{wildcards.dataset}"]
    wildcard_constraints:
        dataset = "energy-balance|hh-end-use"
    conda: "../envs/shell.yaml"
    output: "results/eurostat-{dataset}.tsv.gz"
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


# Swiss data
rule download_ch_energy_data:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: internal["data-sources"][f"swiss-{wildcards.dataset}"]
    output: "results/ch-{dataset}.xlsx"
    conda: "../envs/shell.yaml"
    wildcard_constraints:
        dataset = "energy-balance|industry-energy-balance|end-use"
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


# synthesize eurostat and Swiss data
rule annual_energy_balances:
    message: "Get annual energy balances from Eurostat"
    input:
        energy_balance = "results/eurostat-energy-balance.tsv.gz",
        ch_energy_balance = "results/ch-energy-balance.xlsx",
        ch_industry_energy_balance = "results/ch-industry-energy-balance.xlsx",
        cat_names = workflow.source_path("../resources/energy-balance-category-names.csv"),
        carrier_names = workflow.source_path("../resources/energy-balance-carrier-names.csv")
    output: temp("results/annual-energy-balances.csv")
    params:
        first_year = 2000
    conda: "../envs/default.yaml"
    script: "../scripts/annual_energy_balance.py"


# some source for population. for now, only national resolution
rule download_potentials:
    message: "Download potential data."
    params: url = internal["data-sources"]["potentials"]
    output: temp("results/raw-potentials.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"

rule potentials:
    message: "Unzip potentials."
    input: rules.download_potentials.output[0]
    shadow: "minimal"
    output: "results/population.csv",
    conda: "../envs/shell.yaml"
    shell: "unzip -p {input} national/population.csv > {output}"


rule units_without_shape:
    message: "Dataset of units without geo information."
    input:
        units = "results/{shapes}/shapes.geojson"
    output: "results/{shapes}/units.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/nogeo.py"


# synthesizing eurostat, Swiss data and JRC IDEES
rule annual_heat_demand:
    message: "Calculate national heat demand for household and commercial sectors"
    input:
        hh_end_use = "results/eurostat-hh-end-use.tsv.gz",
        ch_end_use = "results/ch-end-use.xlsx",
        energy_balance = rules.annual_energy_balances.output[0],
        commercial_demand = "results/jrc-idees/tertiary/processed.csv",
        carrier_names = workflow.source_path("../resources/energy-balance-carrier-names.csv")
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
        annual_demand = rules.annual_heat_demand.output["total_demand"],
        electricity = rules.annual_heat_demand.output["electricity"],
        locations = "results/{shapes}/units.csv",
        populations = "results/population.csv"
    conda: "../envs/default.yaml"
    output:
        total_demand = "results/{shapes}/annual-heat-demand-twh.csv",
        electricity = "results/{shapes}/annual-historic-electrified-heat-demand-twh.csv",
    script: "../scripts/rescale.py"
