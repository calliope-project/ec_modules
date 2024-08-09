"""Wrapper functionality for `tsam`'s `TimeSeriesAggregation` methods."""

import pandas as pd
from tsam.timeseriesaggregation import TimeSeriesAggregation

frames = [pd.read_csv(i, index_col=0) for i in snakemake.input]
raw_timeseries = pd.concat(frames, axis=1, verify_integrity=True, copy=False)


aggregation = TimeSeriesAggregation(
    raw_timeseries,
    resolution=snakemake.params.get("resolution"),
    noTypicalPeriods=snakemake.params.get("noTypicalPeriods", 10),
    noSegments=snakemake.params.get("noSegments", 10),
    hoursPerPeriod=snakemake.params.get("hoursPerPeriod", 24),
    clusterMethod=snakemake.params.get("clusterMethod", "hierarchical"),
    evalSumPeriods=snakemake.params.get("evalSumPeriods", False),
    sortValues=snakemake.params.get("sortValues", False),
    sameMean=snakemake.params.get("sameMean", False),
    rescaleClusterPeriods=snakemake.params.get("rescaleClusterPeriods", True),
    weightDict=snakemake.params.get("weightDict"),
    segmentation=snakemake.params.get("segmentation", False),
    extremePeriodMethod=snakemake.params.get("extremePeriodMethod", "None"),
    representationMethod=snakemake.params.get("representationMethod"),
    representationDict=snakemake.params.get("representationDict"),
    distributionPeriodWise=snakemake.params.get("distributionPeriodWise", True),
    segmentRepresentationMethod=snakemake.params.get("segmentRepresentationMethod"),
    predefClusterOrder=snakemake.params.get("predefClusterOrder"),
    predefClusterCenterIndices=snakemake.params.get("predefClusterCenterIndices"),
    solver=snakemake.params.get("solver", "highs"),
    roundOutput=snakemake.params.get("roundOutput"),
    addPeakMin=snakemake.params.get("addPeakMin"),
    addPeakMax=snakemake.params.get("addPeakMax"),
    addMeanMin=snakemake.params.get("addMeanMin"),
    addMeanMax=snakemake.params.get("addMeanMax"),
)

typical_periods = aggregation.createTypicalPeriods()
typical_periods.to_csv(snakemake.output.typical_periods)

# Optional outputs
accuracy_indicators = snakemake.output.get("accuracy_indicators")
if accuracy_indicators:
    aggregation.accuracyIndicators().to_csv(accuracy_indicators)

index_matching = snakemake.output.get("index_matching")
if index_matching:
    aggregation.indexMatching().to_csv(index_matching)
