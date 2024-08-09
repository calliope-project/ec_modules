rule global_regions:
    output:
        shapefile = "resources/global_shapes.geojson",
        plot_shapefile = "results/plots/global_map.png"
    params:
        resolution = "10m",
        category = "cultural",
        name = "admin_1_states_provinces"
    wrapper: "file:../../wrappers/cartopy/shapes/natural_earth"


rule crop_regions:
    input:
        shapefile = rules.global_regions.output.shapefile
    output:
        shapefile = "results/regions/{region}.geojson"
    params:
        column = config["regions"]["column_in_shapefile"],
        name = lambda wc: wc.get("region"),
        bounds = config["regions"]["bounds"]
    conda: "../envs/regions.yaml"
    script: "../scripts/get_regions.py"


