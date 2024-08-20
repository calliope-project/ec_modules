"""Use atlite functions to get renewable profiles."""

rule cutout:
    input:
        shapefile = "results/regions/{region}.geojson"
    output:
        cutout = "results/cutouts/{region}.nc",
        plot_cutout = "results/plots/{region}_cutout.png"
    params:
        time = "2019-05-02",
        features = ["height", "wind", "influx", "temperature", "runoff"],
        offset_degrees = 1,
        module = ["era5"],
    threads: 4
    wrapper: "file:../../wrappers/atlite/cutout-prepare"


rule availability:
    input:
        cutout = "results/cutouts/{region}.nc",
        shapefile = "results/regions/{region}.geojson",
        rasters = ["resources/CORINE2018_V2020_20u1.tif"],
    output:
        availability_matrix = "results/availability/{region}.nc",
        plot_availability_shape = "results/plots/{region}_availability_shape.png",
        plot_availability_matrix = "results/plots/{region}_availability_matrix.png"
    params:
        exclusion_crs = 3035,
        exclusion_resolution = 100,
        shapefile_name_column = config["regions"]["column_in_shapefile"],
        raster_codes = [[12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32]],
        raster_kwargs = [{"invert": False}]
    threads: 4
    wrapper: "file:../../wrappers/atlite/availability"
