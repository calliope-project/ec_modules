# We recommend adding rules that download necessary files here.
if config["use_default_user_resources"]:
    rule download_units:
        message: "Download spatial units."
        params:
            url = internal["resources"]["default_user_shapes"],
        output: "resources/user/shapes_ehighways.geojson"
        conda: "../envs/shell.yaml"
        localrule: True
        shell: "curl -sSLo {output} '{params.url}'"

rule download_potentials:
    message: "Download potential data."
    params: url = internal["resources"]["potentials"]
    output: temp("resources/automatic/raw-potentials.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"

rule unzip_potentials:
    message: "Unzip potentials."
    input: rules.download_potentials.output[0]
    output:
        population = "resources/automatic/{resolution}/population.csv",
        land_cover = "resources/automatic/{resolution}/land-cover.csv"
    conda: "../envs/shell.yaml"
    shadow: "minimal"
    localrule: True
    shell:
        """
        unzip -p {input} '{wildcards.resolution}/population.csv' > '{output.population}'
        unzip -p {input} '{wildcards.resolution}/land-cover.csv' > '{output.land_cover}'
        """

rule download_biofuel_potentials_and_costs:
    message: "Download raw biofuel potential and cost data."
    params: url = internal["resources"]["biofuel_potentials_and_costs"]
    output: "resources/automatic/raw_biofuel_potentials_and_costs.xlsx"
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"
