"""Wrapper for powerplantmatching's powerplant function.

Depending on the given configuration, this may:

- Download the dataset from ppm's repo.
- Trigger the whole powerplantmatching process for:
    - The default ppm configuration.
    - A custom configurations with specific input files.

See ppm's documentation and code for more information.
https://github.com/PyPSA/powerplantmatching
"""

import matplotlib.pyplot as plt
import powerplantmatching as ppm

config = ppm.get_config(
    snakemake.input.get("config"), **snakemake.params.get("config_update", {})
)

plants = ppm.powerplants(
    config=config,
    update=snakemake.params.get("update", False),
    from_url=snakemake.params.get("from_url", False),
    extend_by_vres=snakemake.params.get("extend_by_vres", False),
    extend_by_kwargs=snakemake.params.get("extend_by_kwargs", {}),
    fill_geopositions=snakemake.params.get("fill_geopositions", True),
    filter_missing_geopositions=snakemake.params.get(
        "filter_missing_geopositions", True
    ),
    **snakemake.params.get("collection_kwargs", {}),
)
plants.to_csv(snakemake.output.get("powerplants"))


stats_aggr = snakemake.output.get("stats_aggregated")
if stats_aggr:
    stats = ppm.data.Capacity_stats()
    ppm.plot.fueltype_totals_bar([plants, stats], keys=["Processed", "Statistics"])
    plt.savefig(stats_aggr)
