import matplotlib.pyplot as plt
import pandas as pd


def plot_potentials(path_demand_electricity, path_output):
    """Plot the electricity demand."""
    demand_electricity = pd.read_csv(
        path_demand_electricity, index_col=0, parse_dates=True
    )

    demand_electricity_normed = demand_electricity / demand_electricity.sum(0)

    fig, ax = plt.subplots(figsize=(12, 5), dpi=120)
    demand_electricity_normed.plot(ax=ax, color="blue", alpha=0.01, legend=False)
    ax.set_title(
        "Annual electricity demand,\nnormalized to annual sum for each region."
    )
    plt.savefig(path_output, dpi=120, bbox_inches="tight")


if __name__ == "__main__":
    plot_potentials(snakemake.input[0], snakemake.output[0])
