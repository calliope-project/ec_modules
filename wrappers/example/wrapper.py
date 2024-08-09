"""An example of an append-at-end module for text files."""

import src.appender as appender

if __name__ == "__main__":
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    append_text = snakemake.params.append_text
    capitalize = snakemake.params.capitalize

    assert isinstance(append_text, str)
    assert isinstance(capitalize, bool)
    assert input_file.endswith(".txt")
    assert output_file.endswith(".txt")

    appender.copy_append(input_file, output_file, append_text, capitalize=capitalize)
