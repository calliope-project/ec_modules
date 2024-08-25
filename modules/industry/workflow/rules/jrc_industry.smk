"""JRC-IDEES industry processing rules."""


rule jrc_idees_unzip_industry:
    message: "Unzip {wildcards.country_code} JRC-IDEES industry data."
    input: "resources/jrc_idees/{country_code}.zip",
    output: "resources/jrc_idees/unprocessed_industry_{country_code}.xlsx"
    conda: "../envs/shell.yaml"
    shell: "unzip -p {input} JRC-IDEES-2015_Industry_{wildcards.country_code}.xlsx > {output}"


# TODO: this probably should be done per-country. For missing countries, we could just configure
# the imputation countries in the internal config.
rule jrc_idees_process_industry:
    message: "Process {wildcards.dataset} JRC-IDEES industry data."
    input:
        data = expand(
            "resources/jrc_idees/unprocessed_industry_{country_code}.xlsx",
            country_code=internal["resources"]["jrc_idees"]["spatial_scope"]
        )
    output: "results/processed_industry_{dataset}.nc"
    wildcard_constraints:
        dataset = "energy|production"
    conda: "../envs/default.yaml"
    script: "../scripts/jrc_idees_process_industry.py"
