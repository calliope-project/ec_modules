"""Rules related to wind and solar."""

rule calculate_cf_land:
    message:
        "Generate capacityfactor time series disaggregated by location on "
        "{wildcards.resolution} resolution for {wildcards.technology}."
    input:
        locations="resources/user/shapes_{resolution}.geojson",
        timeseries=ancient("resources/automatic/capacityfactors/{technology}-timeseries.nc"),
        coordinates=ancient("resources/automatic/capacityfactors/wind-onshore-timeseries.nc"),
    params:
        cf_threshold=internal["capacity-factors"]["min"],
        gridcell_overlap_threshold=internal["quality-control"][
            "capacity-factor-gridcell-overlap-threshold"
        ]
    wildcard_constraints:
        technology="|".join(internal["techs"]["land"]),
        year = "|".join([str(i) for i in range(*internal["scope"]["year_range"])])
    output:
        "results/{resolution}/{year}/capacityfactors-{technology}.csv",
    conda:
        "../envs/geo.yaml"
    resources:
        runtime=30,
    script:
        "../scripts/capacityfactors.py"


rule calculate_shared_coast:
    message:
        "Determine share of coast length between EEZ and {wildcards.resolution} units using {threads} threads."
    input:
        units="resources/user/shapes_{resolution}.geojson",
        eez=rules.clip_eez.output[0],
    params:
        polygon_area_share_threshold=internal["quality-control"][
            "shared-coast-polygon-area-share-threshold"
        ]
    output:
        "results/{resolution}/shared-coast.csv",
    threads: 4
    conda:
        "../envs/geo.yaml"
    script:
        "../scripts/shared_coast.py"


rule calculate_cf_offshore:
    message:
        "Generate capacityfactor time series disaggregated by location on "
        "{wildcards.resolution} resolution for wind-offshore."
    input:
        eez=rules.clip_eez.output[0],
        shared_coast="results/{resolution}/shared-coast.csv",
        timeseries=ancient("resources/automatic/capacityfactors/wind-offshore-timeseries.nc"),
    params:
        cf_threshold=internal["capacity-factors"]["min"],
        gridcell_overlap_threshold=internal["quality-control"][
            "capacity-factor-gridcell-overlap-threshold"
        ],
    wildcard_constraints:
        year = "|".join([str(i) for i in range(*internal["scope"]["year_range"])])
    output:
        "results/{resolution}/{year}/capacityfactors-wind-offshore.csv",
    conda:
        "../envs/geo.yaml"
    resources:
        runtime=30,
    script:
        "../scripts/capacityfactors_offshore.py"


rule calculate_area_limits:
    message: "Use technology densities to convert wind & solar {wildcards.resolution} available area to capacity limits."
    input:
        units = "resources/user/shapes_{resolution}.geojson",
        land_eligibility_km2 = "resources/automatic/{resolution}/{scenario}/areas.csv",
    params:
        max_power_densities = config["maximum_installable_mw_per_km2"],
        roof_shares = config["roof_share"],
    output:
        rooftop = "results/{resolution}/{scenario}/rooftop-pv.csv",
        offshore = "results/{resolution}/{scenario}/wind-offshore.csv",
        onshore_and_open_field = "results/{resolution}/{scenario}/open-field-pv-and-wind-onshore.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/capacity_limits.py"
