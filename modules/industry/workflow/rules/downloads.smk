"""Resource download rules."""

rule download_eurostat_energy_balances:
    message: "Download Eurostat energy balance data from euro-calliope datasets."
    params:
        url = internal["resources"]["eurostat_energy_balances_url"]
    conda: "../envs/shell.yaml"
    output:
        balances = "resources/balances/eurostat-energy-balance.tsv.gz"
    localrule: True
    shell: "curl -sSLo {output.balances} {params.url}"


rule download_che_energy_data:
    message: "Get Swiss statistics energy datasets."
    params:
        energy_url = internal["resources"]["che"]["energy_balances_url"],
        industry_url = internal["resources"]["che"]["industry_energy_balances_url"]
    output:
        energy = "resources/balances/che_energy_balances.xlsx",
        industry = "resources/balances/che_industry_energy_balances.xlsx"
    conda: "../envs/shell.yaml"
    localrule: True
    shell:
        """
        curl -sSLo {output.energy} {params.energy_url}
        curl -sSLo {output.industry} {params.industry_url}
        """


rule download_jrc_idees_zipped:
    message: "Download {wildcards.country_code} JRC IDEES zip file."
    params: dataset_url = internal["resources"]["jrc_idees"]["dataset_url"]
    output: temp("resources/jrc_idees/{country_code}.zip")
    conda: "../envs/shell.yaml"
    wildcard_constraints:
        country_code = "|".join(internal["resources"]["jrc_idees"]["spatial_scope"])
    localrule: True
    shell: "curl -f -sSLo {output} '{params.dataset_url}/JRC-IDEES-2015_All_xlsx_{wildcards.country_code}.zip'"
