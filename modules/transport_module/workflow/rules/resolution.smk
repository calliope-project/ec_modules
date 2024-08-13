"""Rules to process units and all the input depending on the resolution."""

rule units:
    message: "Download spatial zones."
    params:
        url = config["module-data-sources"]["spatial-zones"],
    output: "results/automatic/units.geojson"
    conda: "../envs/shell.yaml"
    shell: "curl -sSLo {output} '{params.url}'"


rule units_without_shape:
    message: "Generate a dataset of spatial zones with geoinformation removed."  # FIXME: is there any real reason for this to exist? Why not just use the dataset with geoinformation?
    input:
        units = rules.units.output[0]
    output: "results/automatic/units.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/shapes/nogeo.py"
