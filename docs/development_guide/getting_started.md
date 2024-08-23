# Getting started

We assume you already have `conda` or `mamba` installed in your system.
If you don't, We recommend following `mamba`'s [installation advice](https://github.com/mamba-org/mamba), as it is much faster than regular `conda`.

1. Clone this repository and change into it.

    ```shell
    git clone git@github.com:calliope-project/ec_modules.git
    cd ec_modules
    ```

2. Create a `conda` environment with our setup and activate it.

    ```shell
    conda env create -f environment.yml
    conda activate ec_modules
    ```

3. You are ready to go!
Please look into our [code conventions](conventions.md#code-conventions) and our requirements for developing [modules](modules.md) and [wrappers](wrappers.md) for more details.
