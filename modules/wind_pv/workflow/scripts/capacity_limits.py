import math

import geopandas as gpd
import pandas as pd


def set_capacity_limits(
    path_to_units,
    path_to_land_eligibility_km2,
    roof_shares,
    max_power_densities,
    path_to_onshore_and_open_field_output,
    path_to_offshore_output,
    path_to_rooftop_output,
):
    """Use study data to estimate capacity limit per area."""
    capacities = _from_area_to_installed_capacity(
        land_eligibility_km2=pd.read_csv(path_to_land_eligibility_km2, index_col=0),
        roof_shares=roof_shares,
        maximum_installable_power_density=max_power_densities,
    )
    units = gpd.read_file(path_to_units).set_index("id")
    capacities = capacities.reindex(units.index).fillna(0)
    capacities[
        [
            "eligibility_onshore_wind_monopoly_mw",
            "eligibility_onshore_wind_km2",
            "eligibility_onshore_wind_and_pv_km2",
        ]
    ].to_csv(path_to_onshore_and_open_field_output)
    capacities[
        ["eligibility_offshore_wind_mw", "eligibility_offshore_wind_km2"]
    ].to_csv(path_to_offshore_output)
    capacities[
        [
            "eligibility_rooftop_pv_mw",
            "eligibility_rooftop_pv_n_mw",
            "eligibility_rooftop_pv_e_w_mw",
            "eligibility_rooftop_pv_s_flat_mw",
        ]
    ].to_csv(path_to_rooftop_output)


def _from_area_to_installed_capacity(
    land_eligibility_km2, roof_shares, maximum_installable_power_density: dict
):
    # Roof shares must add to 1
    assert math.isclose(sum(roof_shares.values()), 1, rel_tol=0.01)

    cap = land_eligibility_km2.copy()
    factor_rooftop = maximum_installable_power_density[
        "pv_on_flat_areas"
    ] * roof_shares["flat"] + maximum_installable_power_density[
        "pv_on_tilted_roofs"
    ] * (1 - roof_shares["flat"])
    factor_onshore = maximum_installable_power_density["onshore_wind"]
    factor_offshore = maximum_installable_power_density["offshore_wind"]
    cap["eligibility_rooftop_pv_mw"] = (
        cap["eligibility_rooftop_pv_km2"] * factor_rooftop
    )
    cap["eligibility_offshore_wind_mw"] = (
        cap["eligibility_offshore_wind_km2"] * factor_offshore
    )
    cap["eligibility_onshore_wind_monopoly_mw"] = (
        cap["eligibility_onshore_wind_km2"] * factor_onshore
    )
    cap["eligibility_rooftop_pv_n_mw"] = (
        cap["eligibility_rooftop_pv_km2"]
        * roof_shares["north"]
        * maximum_installable_power_density["pv_on_tilted_roofs"]
    )
    cap["eligibility_rooftop_pv_e_w_mw"] = (
        cap["eligibility_rooftop_pv_km2"]
        * (roof_shares["east"] + roof_shares["west"])
        * maximum_installable_power_density["pv_on_tilted_roofs"]
    )
    cap["eligibility_rooftop_pv_s_flat_mw"] = (
        cap["eligibility_rooftop_pv_km2"]
        * roof_shares["south"]
        * maximum_installable_power_density["pv_on_tilted_roofs"]
    ) + (
        cap["eligibility_rooftop_pv_km2"]
        * roof_shares["flat"]
        * maximum_installable_power_density["pv_on_flat_areas"]
    )
    return cap


if __name__ == "__main__":
    set_capacity_limits(
        path_to_units=snakemake.input.units,
        path_to_land_eligibility_km2=snakemake.input.land_eligibility_km2,
        roof_shares=snakemake.params["roof_shares"],
        max_power_densities=snakemake.params["max_power_densities"],
        path_to_onshore_and_open_field_output=snakemake.output.onshore_and_open_field,
        path_to_offshore_output=snakemake.output.offshore,
        path_to_rooftop_output=snakemake.output.rooftop,
    )
