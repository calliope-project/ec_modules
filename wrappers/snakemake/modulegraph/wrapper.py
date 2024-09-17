"""Wrapper for snakemake - modulegraph."""

__author__ = "Ivan Ruiz Manuel"
__copyright__ = "Copyright 2024, Ivan Ruiz Manuel"
__email__ = "i.ruizmanuel@tudelft.nl"
__license__ = "MIT"

from ec_utils.modules import write_snakemake_modulegraph_png

write_snakemake_modulegraph_png(
    snakemake_dotfile=snakemake.input.dotfile,
    output_path=snakemake.output.pngfile,
    prefixes=snakemake.params.prefixes
)
