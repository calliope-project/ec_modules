"""Snakemake rules to create shapes of administrative regions and exclusive economic zones."""

rule eez:
    message: "Clip exclusive economic zones to study area."
    input: rules.download_eez.output[0]
    output: "build/data/eez.geojson"
    params:
        bounds="{x_min},{y_min},{x_max},{y_max}".format(**config["scope"]["spatial"]["bounds"]),
        countries=",".join(["'{}'".format(country) for country in config["scope"]["spatial"]["countries"]]),
    conda: "../envs/geo.yaml"
    shadow: "minimal"
    shell:
        """
        fio cat --bbox {params.bounds} "zip://{input}"\
        | fio filter "f.properties.TERRITORY1 in [{params.countries}]"\
        | fio collect > {output}
        """
