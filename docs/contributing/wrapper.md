# Wrapper requirements

## File structure

We follow the same file structure guidelines as the [`snakemake-wrappers`](https://snakemake-wrappers.readthedocs.io/en/stable/contributing.html) project.

??? example "Example of a wrapper"
    Most wrappers should follow a structure similar to the one below:

    ```txt
    wrappers/example/
    ├── environment.yaml
    ├── meta.yaml
    ├── README.md
    ├── wrapper.py
    └── test/
        ├── data.csv
        └── Snakefile
    ```

### Obligatory components

- An `environment.yaml` file. A regular `conda` environment with the dependencies of the wrapper. Should be as small as possible.
- A `meta.yaml` file. Documents the `input`, `output` and `params` necessary to use the wrapper.
- A `wrapper.py` file. This is the main wrapper script. It should be concise (no more than ~60 lines of code).
- A `test/` folder. Contains a `Snakefile` with tests and as little data as necessary for executing it.

    !!! warning
        Please avoid placing large files in wrapper tests!

### Optional components

- A `README.md` file with a simple [`mermaid`](https://mermaid.js.org/) diagram showing the wrapper's Input-Output structure.

    ??? example "Wrapper IO diagram"
        ```mermaid
            flowchart LR
            I1(data1.csv) --> W((timeseries))
            I2(data2.csv) --> W
            I3(data3.csv) --> W
            W --> O1(typical_periods.csv)
            W --> |Optional| O2(accuracy_indicators.csv)
            W --> |Optional| O3(index_matching.csv)
        ```
