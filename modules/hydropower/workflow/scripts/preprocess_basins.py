"""Hydrobasin processing and shape 'ribbon' removal."""

import geopandas as gpd

DRIVER = "GPKG"


def preprocess_basins(path_to_basins, path_to_output):
    """Filter and fix basin shapes."""
    basins = gpd.read_file(path_to_basins)
    basins.geometry = basins.geometry.map(_buffer_if_necessary)
    basins.to_file(path_to_output, driver=DRIVER)


# TODO: replace with shapely.make_valid (requires updating many geo.yaml env dependencies to update to shapely 1.8.2)
def _buffer_if_necessary(shape):
    """Fix the basins shapes which are invalid.

    Following the advice given here:
    https://github.com/Toblerity/Shapely/issues/344
    """
    if not shape.is_valid:
        shape = shape.buffer(0.0)
    assert shape.is_valid
    return shape


if __name__ == "__main__":
    preprocess_basins(
        path_to_basins=snakemake.input.basins, path_to_output=snakemake.output[0]
    )
