"""Resource download rules."""


rule download_eurostat:
    message: "Download Eurostat '{wildcards.dataset}' data from euro-calliope datasets."
    params:
        url = lambda wc: internal["resources"]["eurostat"][f"{wc.dataset}"]
    wildcard_constraints:
        country_code = "|".join(internal["resources"]["eurostat"])
    conda: "../envs/shell.yaml"
    output: "resources/eurostat/{dataset}.tsv.gz"
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


rule download_che:
    message: "Download Swiss statistics '{wildcards.dataset}' data."
    params:
        url = lambda wc: internal["resources"]["che"][f"{wc.dataset}"]
    wildcard_constraints:
        country_code = "|".join(internal["resources"]["che"])
    conda: "../envs/shell.yaml"
    output: "resources/che/{dataset}.xlsx"
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


rule download_nuts:
    message: "Download NUTS 2006 shapes."
    params: url = internal["resources"]["eurostat"]["nuts_2006"]
    conda: "../envs/shell.yaml"
    output: "resources/eurostat/nuts_2006.geojson"
    localrule: True
    shell: "curl -sSLo {output} {params.url}"


rule download_jrc_idees_zipped:
    message: "Download {wildcards.country_code} JRC IDEES zip file."
    params: dataset = internal["resources"]["jrc_idees"]["dataset"]
    output: temp("resources/jrc_idees/{country_code}.zip")
    conda: "../envs/shell.yaml"
    wildcard_constraints:
        country_code = "|".join(internal["resources"]["jrc_idees"]["spatial_scope"])
    localrule: True
    shell: "curl -f -sSLo {output} '{params.dataset}/JRC-IDEES-2015_All_xlsx_{wildcards.country_code}.zip'"
