"""Small shapefile creation."""

import geopandas as gpd
from shapely import box

smk = snakemake  # type: ignore


def get_regions(
    shapes: gpd.GeoDataFrame, column: str, name: str, bounds: list
) -> gpd.GeoDataFrame:
    """Extracts regions from a shapefile.

    Args:
        shapes (gpd.GeoDataFrame): Dataframe to cut from.
        column (str): Column to use for filtering.
        name (str | list): Names to select.
        bounds (list): Limit extraction to a zone. Defaults to None.

    Returns:
        gpd.GeoDataFrame: Dataframe with the selected regions.
    """
    assert shapes.crs.is_geographic
    bbox = box(*bounds)
    assert box(*shapes.total_bounds).contains(bbox)
    clipped = shapes.clip(bbox)
    return clipped[clipped[column] == name]


if __name__ == "__main__":
    regions = get_regions(
        shapes=gpd.read_file(smk.input.shapefile),
        column=smk.params.column,
        name=smk.params.name,
        bounds=smk.params.bounds,
    )
    regions.to_file(smk.output.shapefile, driver="GeoJSON")
