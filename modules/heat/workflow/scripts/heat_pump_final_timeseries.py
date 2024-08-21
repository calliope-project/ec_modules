import pandas as pd
import xarray as xr


def _group_end_uses(
    resolution_specific_data: xr.Dataset, annual_demand: xr.DataArray
) -> xr.DataArray:
    weighted_average_da = (
        resolution_specific_data.to_array("end_use")
        .groupby("time.year")
        .apply(_end_use_weighted_ave, annual_demand=annual_demand)
    )
    return weighted_average_da


def prepare_annual_demand(annual_demand: pd.Series) -> xr.DataArray:
    """Restructure annual demand MultiIndex series into a multi-dimensional array.

    Result sums over all building categories and only contains hot water and space
    heating demands (without cooking).
    """
    return (
        annual_demand.rename_axis(columns="id")
        .stack()
        .to_xarray()
        .sel(end_use=["space_heat", "hot_water"])
        .sum("cat_name")
    )


def _end_use_weighted_ave(
    one_year_profile: xr.DataArray, annual_demand: xr.DataArray
) -> xr.DataArray:
    """Calculate weighted average of all heat energy end uses.

    Uses annual demands per spatial unit as the weights.
    """
    year = one_year_profile.time.dt.year[0]
    normalised_demand = annual_demand / annual_demand.sum("end_use")
    demand = one_year_profile * normalised_demand.sel(year=year).drop("year")
    return demand.sum("end_use")


if __name__ == "__main__":
    timeseries_data = xr.open_dataset(snakemake.input.timeseries_data)
    annual_demand = pd.read_csv(snakemake.input.annual_demand, index_col=[0, 1, 2])
    annual_demand_ds = prepare_annual_demand(annual_demand)

    timeseries_data_group_end_use = _group_end_uses(timeseries_data, annual_demand_ds)

    final_df = timeseries_data_group_end_use.astype("float32").to_series().unstack("id")
    final_df.to_csv(snakemake.output[0])
