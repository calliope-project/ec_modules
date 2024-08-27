import pandas as pd
import xarray as xr
from categories.utils import jrc_idees_parser as jrc
from eurocalliopelib import utils as ec_utils


def get_bau_electricity(energy_balances: pd.DataFrame, cat_names: pd.DataFrame) -> xr.DataArray:
    electricity_bau = (
        energy_balances.xs("E7000", level="carrier_code")
        .unstack(["unit", "country", "year"])
        .groupby(cat_names.jrc_idees.dropna().to_dict())
        .sum(min_count=1)
        .stack(["unit", "country"])
        .bfill(axis=1)  # missing old data (BIH and MNE) backfilled
        .stack()
        .apply(ec_utils.tj_to_twh)
        .droplevel(level="unit")
    )

    electricity_bau = pd.concat(
        {"Electricity": electricity_bau}, names=["carrier_name"]
    )
    electricity_bau.index = electricity_bau.index.rename(
        {"cat_code": "cat_name", "country": "country_code"}
    )
    electricity_bau = xr.DataArray.from_series(electricity_bau)
    electricity_bau = jrc.standardize(electricity_bau, "twh", "energy_demand")

    return electricity_bau


def main():
    """Calculate historical electricity demand per country."""
    energy_balances = pd.read_csv(
        snakemake.input.energy_balances,
        index_col=["cat_code", "carrier_code", "unit", "country", "year"],
    ).squeeze("columns")
    cat_names = pd.read_csv(snakemake.input.cat_names, header=0, index_col=0)

    demand = get_bau_electricity(energy_balances, cat_names)
    demand.to_netcdf(snakemake.output[0])



if __name__ == "__main__":
    main()
