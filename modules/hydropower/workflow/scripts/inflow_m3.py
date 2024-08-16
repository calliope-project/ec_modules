"""Use atlite to calculate water inflow timeseries."""

import atlite
import geopandas as gpd
import pandas as pd
import pycountry
import xarray as xr
from shapely.geometry import Point


def determine_water_inflow(
    path_to_cutout,
    path_to_stations,
    path_to_basins,
    first_year,
    final_year,
    path_to_output,
):
    """Determine water inflow timeseries for all plants."""
    plants = _read_plants(path_to_stations)

    inflow_m3 = _water_inflow(plants, path_to_cutout, path_to_basins)
    (
        xr.merge([plants.to_xarray(), inflow_m3])
        .drop("geometry")
        .sel(time=slice(first_year, final_year))
        .to_netcdf(path_to_output)
    )


def _read_plants(path_to_stations):
    plants = pd.read_csv(path_to_stations, index_col=0)
    plants["country_code"] = plants["country_code"].map(
        lambda iso2: pycountry.countries.lookup(iso2).alpha_3
    )
    plants = plants[plants["type"].isin(["HROR", "HDAM"])]
    return gpd.GeoDataFrame(
        plants, geometry=list(map(Point, zip(plants.lon, plants.lat)))
    )


def _water_inflow(plants, path_to_cutout, path_to_basins):
    cutout = atlite.Cutout(path=path_to_cutout)
    inflow = cutout.hydro(plants, path_to_basins).rename(plant="id").rename("inflow_m3")
    return inflow


if __name__ == "__main__":
    determine_water_inflow(
        path_to_cutout=snakemake.input.runoff,
        path_to_stations=snakemake.input.stations,
        path_to_basins=snakemake.input.basins,
        first_year=str(snakemake.wildcards.first_year),
        final_year=str(snakemake.wildcards.final_year),
        path_to_output=snakemake.output[0],
    )
