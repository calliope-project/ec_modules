"""Wrapper for atlite cutout generation."""

import atlite
import xarray as xr

CUTOUT_OPTIONAL = ["dx", "dy", "dt", "chunks"]
PREPARE_OPTIONAL = ["features", "tmpdir", "overwrite", "compression"]

# Customise cutout
cutout_params = {}
input_data = snakemake.input.get("data")
if input_data:
    cutout_params["data"] = xr.open_dataset(input_data)

for param in CUTOUT_OPTIONAL:
    value = snakemake.params.get(param)
    if value is not None:
        cutout_params[param] = value

offset = snakemake.params.get("cutout_offset", 0)
x = slice(snakemake.params.x[0] - offset, snakemake.params.x[1] + offset)
y = slice(snakemake.params.y[0] - offset, snakemake.params.y[1] + offset)

cutout = atlite.Cutout(
    path=snakemake.output.path,
    x=x,
    y=y,
    time=snakemake.params.time,
    module="era5",
    **cutout_params,
)

# Download and prepare data.
prepare_params = {}
for param in PREPARE_OPTIONAL:
    value = snakemake.params.get(param)
    if value is not None:
        prepare_params[param] = value
cutout.prepare(**prepare_params)
