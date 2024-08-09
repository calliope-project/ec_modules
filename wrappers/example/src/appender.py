"""Appender functionality, independent of Snakemake's mess!"""

from pathlib import Path
from shutil import copy


def copy_append(
    input_path: str, output_path: str, text: str, capitalize: bool = False
) -> None:
    """Copy a file and append text to it."""
    from_file = Path(input_path)
    to_file = Path(output_path)

    to_file.unlink(missing_ok=True)
    to_file.parent.mkdir(parents=True, exist_ok=True)
    copy(from_file, to_file)

    with to_file.open(mode="a") as file:
        if capitalize:
            file.write(text.capitalize())
        else:
            file.write(text)
