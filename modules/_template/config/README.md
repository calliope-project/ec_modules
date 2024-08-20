# Module configuration

User-modifiable module settings should go in `default.yaml`.

You should also provide a brief explanation of each setting below.

## Example

- shapefile: `.geojson` file containing at least a `geometry` column, a column with the country code in ISO 3166-1 alpha 3 standard, and a column with unique subregion IDs.
  - path: URL to the shape, local path, or relative path (must be relative to the main workflow).
  - download: if `True`, the module will use the `path` to download the shapefile.
  - country-id-column: name of the column containing ISO country codes in alpha-3 format.
  - subregion-id-column: name of the column containing unique subregion identifiers.
- year-slice: years to process, in the form [start_year, end_year].  Inclusive (e.g., [2015, 2017] will include 2015, 2016, 2017).
