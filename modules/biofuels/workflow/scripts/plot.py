import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd


def plot_potentials(path_shapes, path_potentials, path_costs, path_output):
    """Plot the biofuel potentials on a map."""
    shapes = gpd.read_file(path_shapes).to_crs(epsg=3035)
    potentials = pd.read_csv(path_potentials)
    with open(path_costs) as file:
        costs = float(file.read())
    potentials = shapes.merge(potentials, on="id")

    fig, ax = plt.subplots(figsize=(5, 5), dpi=120)
    potentials.plot(
        ax=ax,
        column="biofuel_potential_mwh_per_year",
        cmap="Greens",
        edgecolor="black",
        legend=True,
    )
    ax.set_title(
        "Biofuel potential [MWh/year],\n at a cost of {costs:.2f} €/MWh".format(
            costs=costs
        )
    )
    ax.set_axis_off()
    plt.savefig(path_output, dpi=120, bbox_inches="tight")


if __name__ == "__main__":
    plot_potentials(
        snakemake.input.shapes,
        snakemake.input.potentials,
        snakemake.input.costs,
        snakemake.output[0],
    )
