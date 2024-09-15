"""Generate offshore capacityfactor time series."""

import geopandas as gpd
import pandas as pd
import xarray as xr
from eurocalliopelib.geo import (
    area_weighted_time_series,
    convert_old_style_capacity_factor_time_series,
)

EPSG3035 = "EPSG:3035"


def capacityfactors(
    path_to_eez,
    path_to_shared_coast,
    path_to_timeseries,
    cf_threshold,
    path_to_result,
    gridcell_overlap_threshold,
    year,
):
    """Generate offshore capacityfactor time series for each location."""
    eez = gpd.read_file(path_to_eez).set_index("MRGID").to_crs(EPSG3035).geometry
    shared_coast = pd.read_csv(path_to_shared_coast, index_col=[0, 1], squeeze=True)

    ts = xr.open_dataset(path_to_timeseries)
    ts = ts.sel(time=slice(year, year))
    # xarray will silently miss the fact that data doesn't exist with slice
    if year not in ts.time.to_index().year.astype(str).unique():
        raise ValueError(
            f"Cannot access capacity factor data for timeseries {path_to_timeseries} "
            f"for {year}."
        )
    ts = convert_old_style_capacity_factor_time_series(ts)

    capacityfactors_per_eez = area_weighted_time_series(
        shapes=eez,
        spatiotemporal=ts,
        gridcell_overlap_threshold=gridcell_overlap_threshold,
    )
    capacityfactors = _allocate_to_onshore_locations(
        capacityfactors_per_eez, shared_coast
    )
    capacityfactors.where(capacityfactors >= cf_threshold, 0).to_csv(path_to_result)


def _allocate_to_onshore_locations(capacityfactors_per_eez, shared_coast):
    shared_coast_norm_per_id = shared_coast.groupby("id").apply(lambda x: x / x.sum())
    return (
        capacityfactors_per_eez.rename_axis(index={"timestep": "timesteps"}, columns="MRGID")
        .mul(shared_coast_norm_per_id, axis=1)
        .groupby("id", axis=1)
        .sum()
    )


if __name__ == "__main__":
    capacityfactors(
        path_to_eez=snakemake.input.eez,
        path_to_shared_coast=snakemake.input.shared_coast,
        path_to_timeseries=snakemake.input.timeseries,
        cf_threshold=float(snakemake.params.cf_threshold),
        gridcell_overlap_threshold=float(snakemake.params.gridcell_overlap_threshold),
        year=str(snakemake.wildcards.year),
        path_to_result=snakemake.output[0],
    )
