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

3. Create your module or wrapper using one of our standard [`copier`](https://github.com/copier-org/copier) templates.

    ??? example "Example: using a module template"

        With the `ec_modules` environment activated, type:

        ```shell
        copier copy modules/_template/ modules/
        ```

        You'll be promted with some questions. After answering them, `copier` will auto-generate the module for you!

        ```html
        🎤 What is your module's name?
        wind_offshore
        🎤 Please give a brief sentence describing your module.
        A module to estimate offshore-wind potentials for arbitrary subregions in Europe.
        🎤 We auto-generate an MIT license for you. Please provide your full name.
        E. Dantès
        🎤 We auto-generate an MIT license for you. Please provide the name of your institution.
        Morrel Technical Institute
        🎤 We auto-generate an MIT license for you. Please provide an email address.
        e.dantes@morrel-ti.edu
        🎤 We auto-generate an MIT license for you. What year is this?
        1815
        ```

    ??? example "Example: using a wrapper template"

        Similar to the module example, call `copier` with the following:

        ```shell
        copier copy wrappers/_template/ wrappers/
        ```

        In the case of wrappers, there are some additional answers you must provide.

        ```html
        🎤 What is the name of the tool you are designing a wrapper for?
        gregor
        🎤 Please provide a valid link to the tool's official website.
        https://github.com/jnnr/gregor
        🎤 What is the name of the wrapper?
        snip
        🎤 Please give a brief sentence describing your wrapper.
        Snip a raster file into a smaller raster file.
        🎤 We auto-generate an MIT license for you. Please provide your full name.
        G. Samsa
        🎤 We auto-generate an MIT license for you. Please provide the name of your institution.
        Bekannt University
        🎤 We auto-generate an MIT license for you. Please provide an email address.
        g.s@un.bekannt.edu
        🎤 We auto-generate an MIT license for you. What year is this?
        1915
        ```

You are ready to go!
Please look into our [code conventions](conventions.md#code-conventions) and our requirements for developing [modules](modules.md) and [wrappers](wrappers.md) for more details.
