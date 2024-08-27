"""Process Chemical industry according to the method in Sector-Coupled Euro-Calliope."""

from typing import TypedDict

import pandas as pd
import xarray as xr
from utils import filling
from utils import jrc_idees_parser as jrc

H2_LHV_KTOE = 2.863  # 0.0333 TWh/kt LHV -> 2.863ktoe/kt
MWH_PER_T_TO_KTOE_PER_KT = 0.08598  # 1 MWh/t -> 1 GWh/kt -> 0.08598 ktoe/kt
METHANOL_LHV_KTOE = 0.476  # 19.915 MJ/kg LHV -> 19915000MJ/kt -> 0.476ktoe/kt
KTOE_TO_TWH = 0.01163
KT_TO_T = 1000

CAT_NAME_CHEMICALS = "Chemicals Industry"

Demand = TypedDict("Demand", {"energy": xr.DataArray, "mass": xr.DataArray})


def get_chemicals_demand_df(
    # config: dict,
    energy_balances: str,
    cat_names: str,
    carrier_names: str,
    jrc_industry_energy: str,
    jrc_industry_production: str
) -> Demand:
    """Execute the data processing pipeline for the "Chemicals Industry" sub-sector.

    Args:
        config (dict): chemicals industry sector configuration.
        energy_balances (str): country energy balances (usually from eurostat).
        cat_names (str): eurostat category mapping file.
        carrier_names (str): eurostat carrier name mapping file.
        jrc_industry_energy (str): jrc country-specific industrial energy demand file.
        jrc_industry_production (str): jrc country-specific industrial production file.
        output_path (str): location of chemicals demand output file. Defaults to None.
    """
    # Load data
    energy_balances_df = pd.read_csv(
        energy_balances,
        index_col=["cat_code", "carrier_code", "unit", "country", "year"],
    ).squeeze("columns")
    cat_names_df = pd.read_csv(cat_names, header=0, index_col=0)
    carrier_names_df = pd.read_csv(carrier_names, header=0, index_col=0)
    jrc_energy = xr.open_dataset(jrc_industry_energy)
    jrc_prod = xr.open_dataarray(jrc_industry_production)
    jrc.check_units(jrc_energy, jrc_prod)

    # Ensure dataframes only have data specific to this industry
    cat_names_df = cat_names_df[cat_names_df["jrc_idees"] == CAT_NAME_CHEMICALS]
    jrc_energy = jrc_energy.sel(cat_name=CAT_NAME_CHEMICALS)
    jrc_prod = jrc_prod.sel(cat_name=CAT_NAME_CHEMICALS)

    # Process data
    new_demand = transform_jrc_subsector_demand(
        jrc_energy,
        jrc_prod,
        # config
    )

    for var in new_demand:
        new_demand[var] = filling.fill_missing_countries_years(
            energy_balances_df, cat_names_df, carrier_names_df, new_demand[var]
        )

    return new_demand


def transform_jrc_subsector_demand(
    jrc_energy: xr.Dataset,
    jrc_prod: xr.Dataset,
    # config_chemicals: dict
) -> Demand:
    # Get the electricity consumption (assuming every demand that could be met by
    # electricity is met by electricity)
    electrical_consumption = (
        (jrc.replace_carrier_final_demand("Electricity", jrc_energy))
        .drop_sel(subsection="Low enthalpy heat")
        .sum(dim=["section", "subsection"])
        .drop_vars("cat_name")
    )

    # Gather relevant params for Basic chemicals in the Chemicals Industry
    h2_demand = {  # t/t
        "Ethylene": 0,
        "Propylene": 0,
        "BTX": 0,
        "Ammonia": 0.178,
        "Methanol": 0.189,
    }
    methanol_demand = {
        "Ethylene": 2.83,
        "Propylene": 2.83,
        "BTX": 4.3,
        "Ammonia": 0,
        "Methanol": 1,
    }
    co2_demand = {  # tCO2/t
        "Ethylene": 0,
        "Propylene": 0,
        "BTX": 0,
        "Ammonia": 0.112,  # inc. 0.35t_urea/t_ammonia
        "Methanol": 1.373,
    }

    energy_demand = {  # MWh/t
        "Ethylene": 1.4,  # MTO process
        "Propylene": 1.4,  # MTO process
        "BTX": 1.4,  # MTA process
        "Ammonia": 2.05,  # inc. 0.35t_urea/t_ammonia
        "Methanol": 1.5,  # auxiliary demand, could just be assumed as already included in existing electricity demand
    }

    mass = {  # Bazzanella and Ausfelder, 2017
        "Ethylene": 21.7,
        "Propylene": 17,
        "BTX": 15.7,
        "Ammonia": 17,
        "Methanol": 2,
    }
    molar_mass = {
        "Ethylene": 28.05,
        "Propylene": 42.08,
        "BTX": 93,
        "Ammonia": 17.01,
        "Methanol": 32.04,
    }
    moles = {k: mass[k] / molar_mass[k] for k in mass}
    molar_ratio = {k: v / sum(moles.values()) for k, v in moles.items()}

    basic_chemicals_mass = jrc_prod.sel(
        produced_material="Basic chemicals (ethylene eq.)"
    )
    basic_chemicals_moles = basic_chemicals_mass / molar_mass["Ethylene"]
    masses = {
        k: basic_chemicals_moles * molar_ratio[k] * molar_mass[k] for k in molar_ratio
    }

    chem_h2_consumption = (
        sum(masses[k] * h2_demand[k] * H2_LHV_KTOE for k in h2_demand) * KTOE_TO_TWH
    )
    chem_h2_consumption = chem_h2_consumption.assign_attrs(units="twh").assign_coords(
        carrier_name="Hydrogen"
    )

    chem_methanol_consumption = (
        sum(masses[k] * methanol_demand[k] * METHANOL_LHV_KTOE for k in methanol_demand)
        * KTOE_TO_TWH
    )
    chem_methanol_consumption = chem_methanol_consumption.assign_attrs(
        units="twh"
    ).assign_coords(carrier_name="Methanol")

    chem_co2_consumption = sum(masses[k] * co2_demand[k] for k in co2_demand) * KT_TO_T
    chem_co2_consumption = chem_co2_consumption.assign_attrs(
        units="tonne"
    ).assign_coords(carrier_name="Carbon Dioxide")

    chem_energy_consumption = (
        sum(
            masses[k] * energy_demand[k] * MWH_PER_T_TO_KTOE_PER_KT
            for k in energy_demand
        )
        * KTOE_TO_TWH
        + electrical_consumption
    )
    chem_energy_consumption = chem_energy_consumption.assign_attrs(
        units="twh"
    ).assign_coords(carrier_name="Electricity")

    # Space heat
    space_heat_demand = (
        jrc_energy["useful"]
        .sel(subsection="Low enthalpy heat")
        .sum(dim=["section", "carrier_name"])
        .drop_vars(["subsection", "cat_name"])
        .assign_attrs(units="twh")
        .assign_coords(carrier_name="Low enthalpy heat")
    )

    # Combine and transform to energy demand
    energy_demand = xr.concat(
        [
            chem_energy_consumption,
            chem_h2_consumption,
            chem_methanol_consumption,
            space_heat_demand,
        ],
        dim="carrier_name",
        coords="minimal",
    )

    # Prettify
    demand = Demand(
        energy=jrc.standardize(energy_demand, "twh", "energy_demand"),
        mass=jrc.standardize(chem_co2_consumption, "tonnes", "co2_demand")
    )

    return demand


if __name__ == "__main__":
    chem_demand = get_chemicals_demand_df(
        # config=snakemake.params.config,
        energy_balances=snakemake.input.energy_balances,
        cat_names=snakemake.input.cat_names,
        carrier_names=snakemake.input.carrier_names,
        jrc_industry_energy=snakemake.input.jrc_industry_energy,
        jrc_industry_production=snakemake.input.jrc_industry_production,
    )

    chem_demand["energy"].to_netcdf(snakemake.output.energy)
    chem_demand["mass"].to_netcdf(snakemake.output.mass)
