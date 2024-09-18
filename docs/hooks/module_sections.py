"""Generate documentation for all the existing modules."""

import subprocess
import tempfile
from pathlib import Path

DOCS_PATH = Path("docs/")
MODULES_PATH = Path("modules/")
TEMPDIR = tempfile.TemporaryDirectory(dir=DOCS_PATH / "modules")

TEMPLATE = """
--8<-- "{readme}:intro"

???+ info "General summary"

    --8<-- "{readme}:mermaid"

??? info "`snakemake` execution steps"

    ![rulegraph]({rulegraph})

## How to use

All you need is three easy steps:

1. Configuring the module correctly.
2. Placing files in the `resources/user` folder with the correct filename (if needed).

    ??? example "User file with wildcard"

        Imagine that a module requests `resources/user/something_{{wildcard}}.txt`,
        with no limitations in wildcard naming.
        The following would be valid:

        ```txt
        resources/
            └── user/
                └── something_foobar.txt
        ```

        For prefixed modules (i.e., using `prefix: some_prefix` in `snakemake`):

        ```txt
        some_prefix/
        └── resources/
            └── user/
                └── something_foobar.txt
        ```

3. Requesting results using the correct filenames and wildcards (if needed):

    ```bash
    snakemake results/{{wildcard1}}/filename_{{wildcard2}}.csv --cores 2
    ```

### Configuration

We recommend to start with the following configuration,
and then experiment with other settings.

???+ example "Default configuration"

    ```yaml
    --8<-- "{default_config}"
    ```

Here you can find a complete description of the module's configuration schema.

??? info "Configuration schema"

    ```yaml
    --8<-- "{schema}"
    ```

### Wildcards

These `snakemake` wildcards should be specified for particular files in
the `resource/user/` and `results/` folders.
They can be used to dynamically alter the behaviour of the module.

--8<-- "{readme}:wildcards"

### User files

These files must be provided to the module beforehand in `resources/user/`.

--8<-- "{readme}:user"

### Result files

These files are the final results of this module.
They will be placed in the `results/` folder when requested by your `snakemake` project.

--8<-- "{readme}:results"

## Attribution

!!! quote "Citation"

    {citation}

!!! quote "References"

    --8<-- "{readme}:references"

??? info "Contributors"

    --8<-- "{authors}:authors"

??? info "License"

    --8<-- "{license}"

"""


def on_pre_build(config):
    """Generate standard documentation page per module.

    Automatically called by mkdocs if the hook is configured.
    """
    for module_dir in MODULES_PATH.iterdir():
        name = module_dir.name
        create_module_docfile(name)


def create_module_docfile(name: str, docs_dir: Path, module_dir: Path):
    """Save a fully documented page for the requested module."""
    authors = module_dir / "AUTHORS.md"
    readme = module_dir / "README.md"
    license = module_dir / "LICENSE.txt"
    default_config = module_dir / "config/default.yaml"
    schema = module_dir / "workflow/schemas/config.schema.yaml"
    rulegraph = create_module_rulegraph(name, Path(TEMPDIR.name), module_dir)
    assert all(
        [
            file.exists()
            for file in [authors, readme, default_config, rulegraph, schema, license]
        ]
    )

    citation = get_apa_citation(module_dir)

    text = TEMPLATE.format(
        readme=str(readme),
        default_config=str(default_config),
        rulegraph=str(rulegraph.relative_to(docs_dir.resolve())),
        citation=citation,
        schema=schema,
        authors=authors,
        license=license,
    )

    with open(docs_dir / f"{name}.md", "w") as file:
        file.write(text)


def create_module_rulegraph(name, prefix: Path, module_dir: Path) -> Path:
    """Run snakemake rulegraph commands and save a .png in the requested location."""
    tmp_rulegraph = prefix / f"{name}.png"
    command = f"snakemake --rulegraph | dot -Tpng > {str(tmp_rulegraph.resolve())}"
    subprocess.run(command, shell=True, check=True, cwd=str(module_dir.resolve()))

    return tmp_rulegraph


def get_apa_citation(module_dir: Path):
    """Return APA citation if the .cff file is correct."""
    dir = str(module_dir)
    subprocess.run("cffconvert --validate", shell=True, check=True, cwd=dir)
    citation = subprocess.run(
        "cffconvert -f apalike", shell=True, check=True, cwd=dir, capture_output=True
    )
    return citation.stdout.decode()


if __name__ == "__main__":
    create_module_docfile(
        "wind_pv", Path("docs/modules/"), Path("modules/wind_pv/")
    )
    pass
