# Wrapper for `atlite` cutout functionality

Downloads and stores a cutout of the requested latitude/longitude and saves it locally for reuse.

>[!important]
>`atlite` data reading is not combinatorial. Three cases are possible AND EXCLUSIVE EACH OTHER:
>
>1. Cutout exists in given path (output.path) -> take it.
>2. Otherwise, input.data (optional) given -> take it.
>3. Otherwise, build cutout using snakemake.params.

## Example

```snakemake
rule atlite_cutout:
    output:
        path = "output/test.nc",
    params:
        time = "2012-01-01",
        x = [3.9, 4.7],
        y = [50, 52.2],
        features = ["wind", "runoff"]
    threads: 4
    wrapper: github("calliope-project/ec_modules", path="wrappers/atlite/cutout")

```
