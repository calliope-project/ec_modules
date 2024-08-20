"""Rules to process units and all the input depending on the resolution."""

# FIXME: is there any real reason for this to exist? Why not just use the dataset with geoinformation?
rule units_without_shape:
    message: "Generate a dataset of spatial zones with geoinformation removed."
    input:
        units = rules.download_units.output[0]
    output: "results/downloads/units.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/shapes/nogeo.py"
