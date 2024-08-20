
if config["shapefile"]["download"]:
    rule download_shapefile:
        message: "Downloading the configured shapefile."
        params:
            url = config["shapefile"]["path"],
        output: "results/downloads/shapes.geojson"
        conda: "../envs/shell.yaml"
        shell: "curl -sSLo {output} '{params.url}'"


rule hello_world:
    message: "I am a module and that's OK!"
    output: "results/hello.txt"
    conda: "../envs/shell.yaml"
    script: "../scripts/example.py"
