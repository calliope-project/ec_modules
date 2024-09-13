import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd


def sort_columns_individually(df):
    """Sort the columns of a DataFrame individually."""
    df = df.copy()
    for col in df:
        df[col] = df[col].sort_values(ascending=False, ignore_index=True)
    return df


def plot_timeseries(path_demand_electricity, path_output):
    """Plot the electricity demand."""
    demand_electricity = pd.read_csv(
        path_demand_electricity, index_col=0, parse_dates=True
    )

    demand_electricity_normed = demand_electricity / demand_electricity.sum(0)
    demand_electricity_normed_sorted = sort_columns_individually(
        demand_electricity_normed.reset_index(drop=True)
    )

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 5), dpi=120)

    demand_electricity_normed.plot(ax=ax1, color="blue", alpha=0.01, legend=False)
    demand_electricity_normed_sorted.plot(ax=ax2, color="red", alpha=0.1, legend=False)
    fig.suptitle(
        "Annual electricity demand,\nnormalized to annual sum for each region."
    )
    plt.savefig(path_output, dpi=120, bbox_inches="tight")


def plot_map(path_demand_electricity, path_shapes, path_output):
    """Plot the electricity demand on a map."""
    demand_electricity = pd.read_csv(
        path_demand_electricity, index_col=0, parse_dates=True
    )
    shapes = gpd.read_file(path_shapes).set_index("id")

    demand_electricity = demand_electricity.sum(0)
    demand_electricity.name = "demand_electricity"
    demand_electricity_shapes = shapes.join(demand_electricity)
    demand_electricity_shapes = demand_electricity_shapes.to_crs(epsg=3035)

    fig, ax = plt.subplots(figsize=(12, 5), dpi=120)
    demand_electricity_shapes.plot(
        ax=ax,
        column="demand_electricity",
        cmap="Blues",
        legend=True,
        edgecolor="black",
        legend_kwds={"label": "Electricity demand"},
    )
    ax.set_axis_off()
    ax.set_title("Annual electricity demand.")
    plt.savefig(path_output, dpi=120, bbox_inches="tight")


if __name__ == "__main__":
    plot_timeseries(snakemake.input.demand_electricity, snakemake.output.timeseries)
    plot_map(
        snakemake.input.demand_electricity,
        snakemake.input.shapes,
        snakemake.output.maps,
    )
