# README
**NOTE:** There exists a newer and more complete version of these notes.
Ask @amotta if you need this information (tag:
2018-08-01-Reconstruction-Workflow-and-State-Files)

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
   See `+connectEM/getAggloQueryOverlapA.m`
2. Mapping of tracings onto axon agglomerates  
   See `+connectEM/getAggloQueryOverlapB.m`
3. Mapping of tracings onto endings of axon agglomerates  
   See `+connectEM/flightEndingOverlapRun.m` and
   `+connectEM/flightEndingOverlap.m`
4. Case distinction  
   See `+connectEM/makeEndingCaseDistinctions.m`
5. Patching in of flight queries
   See `+connectEM/createNewSuperagglos.m`

### Case distinctions

 1. Flight does not reach any axon agglomerate
 2. Flight overlaps with axon agglomerate at end, but not at start
 3. Flight overlaps with axon agglomerate at start, but not at end
    (after exclusion of type 10 flights)
 4. Neither end of flight overlaps with an axon ending
 5. Start, but not end of flight does overlap with an axon ending (after
    exclusion of type 11 flights)
 6. End, but not start of flight does overlap with an axon ending
 7. Start and end of flight overlap with same axon ending
 8. Start and end of flight overlap with same axon agglomerate, but not
    the same ending (i.e., different endings or no endings)
 9. "Correct" flight that overlaps with endings of two different axon
    agglomerates (after removal of type 12 flights)
10. Special case of type 3: Flight overlaps with axon agglomerate at
    start and comes closer than 3 µm to the end of dataset
11. Special case of type 5: Axon ending at start was queried multiple
    times and was answered with type 5 flights that reach different axon
    agglomerates.
12. Special case of type 9: At least one of the involved axon endings
    was reached by case 9 flights that started at different axon
    endings.
13. 
