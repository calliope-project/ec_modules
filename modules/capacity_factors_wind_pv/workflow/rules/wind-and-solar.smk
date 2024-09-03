"""Rules related to wind and solar."""


rule capacity_factors_onshore_wind_and_solar:
    message:
        "Generate capacityfactor time series disaggregated by location on "
        "{wildcards.resolution} resolution for {wildcards.technology}."
    input:
        locations=rules.download_units.output[0],
        timeseries=ancient("resources/capacityfactors/{technology}-timeseries.nc"),
        coordinates=ancient("resources/capacityfactors/wind-onshore-timeseries.nc"),
    params:
        cf_threshold=internal["capacity-factors"]["min"],
        gridcell_overlap_threshold=internal["quality-control"][
            "capacity-factor-gridcell-overlap-threshold"
        ],
        first_year=config["scope"]["temporal"]["first-year"],
        final_year=config["scope"]["temporal"]["final-year"],
        trim_ts=internal["capacity-factors"]["trim-ninja-timeseries"],
    wildcard_constraints:
        technology="wind-onshore|rooftop-pv|open-field-pv|rooftop-pv-n|rooftop-pv-e-w|rooftop-pv-s-flat",
    output:
        "results/{resolution}/capacityfactors-{technology}.csv",
    conda:
        "../envs/geo.yaml"
    resources:
        runtime=30,
    script:
        "../scripts/wind-and-solar/capacityfactors.py"


rule shared_coast:
    message:
        "Determine share of coast length between EEZ and {wildcards.resolution} units using {threads} threads."
    input:
        units=rules.download_units.output[0],
        eez=rules.eez.output[0],
    params:
        polygon_area_share_threshold=internal["quality-control"][
            "shared-coast-polygon-area-share-threshold"
        ],
    output:
        "results/{resolution}/shared-coast.csv",
    threads: 4
    conda:
        "../envs/geo.yaml"
    script:
        "../scripts/wind-and-solar/shared_coast.py"


rule capacity_factors_offshore:
    message:
        "Generate capacityfactor time series disaggregated by location on "
        "{wildcards.resolution} resolution for wind-offshore."
    input:
        eez=rules.eez.output[0],
        shared_coast=rules.shared_coast.output[0],
        timeseries=ancient("resources/capacityfactors/wind-offshore-timeseries.nc"),
    params:
        cf_threshold=internal["capacity-factors"]["min"],
        gridcell_overlap_threshold=internal["quality-control"][
            "capacity-factor-gridcell-overlap-threshold"
        ],
        first_year=config["scope"]["temporal"]["first-year"],
        final_year=config["scope"]["temporal"]["final-year"],
        trim_ts=internal["capacity-factors"]["trim-ninja-timeseries"],
    output:
        "results/{resolution}/capacityfactors-wind-offshore.csv",
    conda:
        "../envs/geo.yaml"
    resources:
        runtime=30,
    script:
        "../scripts/wind-and-solar/capacityfactors_offshore.py"
