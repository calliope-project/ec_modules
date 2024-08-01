"""Wrapper for atlite cutout generation."""

import atlite

OPTIONAL_INPUTS = ["data", "gebco_path"]

offset = snakemake.params.cutout_offset
x = slice(snakemake.params.x[0] - offset, snakemake.params.x[1] + offset)
y = slice(snakemake.params.y[0] - offset, snakemake.params.y[1] + offset)

cutout_kwargs = snakemake.params.get("cutout_kwargs", {})

for optional_input in OPTIONAL_INPUTS:
    value = snakemake.input.get(optional_input)
    if value:
        cutout_kwargs[optional_input] = value

breakpoint()
cutout = atlite.Cutout(
    path=snakemake.output.cutout,
    module=snakemake.params.module,
    x=x,
    y=y,
    time=snakemake.params.time,
    **cutout_kwargs,
)

cutout.prepare(snakemake.params.features, **snakemake.params.get("prepare_kwargs", {}))
