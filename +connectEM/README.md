# README

## Blood vessel detection and correction for intensity gradient
See `preprocessing/preprocessData.m`

## Nuclei and myelin detection
See
* `preprocessing/additionalHeuristics.m`,
* `preprocessing/localDetectionNuclei.m`, and
* `preprocessing/localDetectionMyelin.m`.

## Pipeline

See `pipeline/runPipeline.m`

### Volume segmentation
See
* `bigFwdPass.m`,
* `+Myelin/runFix.m` and `+Myelin/enforceMyelinSegments.m`, and
* `miniSegmentation.m`, `segmentation/segmentForPipeline.m`, and `segmentation/watershedSeg_v1_cortex.m`.

### Segment neighborhood graph
See `graphConstruction.m` and `+SynEM/+Svg/findEdgesAndBorders.m`.

### SynEM classifier
See `+SynEM/+Seg/pipelineRun.m`

### ConnectEM classifier
See
* `+connectEM/calculateFeatures.m`, and
* `+connectEM/predictDataset.m`.

### TypeEM classifier
See
* `+TypeEM/+Pipeline/buildFeatures.m`, and
* `+TypeEM/+Pipeline/buildPredictions.m`

## Automated axon reconstruction
1. Grid search (`+connectEM/aggloGridSearch.m`). Yields
   `aggloGridSearch/search05_00564.mat` based on the following set of
   parameters
     * min. border size: 60 vx
     * min. segment size: 300 vx
     * min. axon type probability: 50 %
     * min. neurite continuity probability: 97 %
2. Grid search (`+connectEM/aggloGridSearchDir.m`). Yields
   `aggloGridSearch6/6_01_00046/metricsFinal.mat` based on the following
   set of parameters
     * min. latent score: 0.8
     * min. border size: 40 vx
     * min. directionality score: 0.9
     * min. neurite continuity probability: 80 %
     * min. axon type probability: 30 %
     * recursion steps: 10
3. Conversion from segment equivalence classes to super-agglomerate
   representation (`+connectEM/aggloPreprocessing.m`). Yields
   `axons_01.mat`.
4. Merging of super-agglomerates that are linked by non-myelin
   cross-segmentation cube edges (`+connectEM/aggloPreprocessing.m`).
   Yields `axons_02.mat`.
5. Split super-agglomerates in proximity to dataset boundary, where
   the decreasing alignment quality results in a heightened rate of
   merge errors (`+connectEM/splitBordersMeta.m`). Yields
   `axons_04.mat`.

## Automated dendrite reconstruction
Grid search (`+connectEM/aggloGridSearch.m`) from commit
618b07f4eaaf5a98cb00d2926dcc4354355c1a7d. Yields
`aggloGridSearch/search03_00514.mat` based on the following set of
parameters:

* min. border size: 300 vx
* min. segment size: 500 vx
* min. dendrite type probability: 30 %
* min. neurite continuity probability: 98 %

## FocusEM ending queries
1. Parsing of NML files and segment overlap  
   See TODO
2. Mapping of tracings onto axon agglomerates  
   See `+connectEM/getAggloQueryOverlapB.m`
3. Mapping of tracings onto endings of axon agglomerates  
   See `+connectEM/flightEndingOverlapRun.m` and
   `+connectEM/flightEndingOverlap.m`
4. Case distinction  
   See `+connectEM/makeEndingCaseDistinctions.m`
5. Patching in of flight queries
   See `+connectEM/createNewSuperagglos.m`
