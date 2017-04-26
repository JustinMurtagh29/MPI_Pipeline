% allows for setting a specific starting point in the pipeline
% if for example some steps were already done before like in the case of an error mid-pipeline 
classdef PipelineStep < uint8
   % step name corresponding to job name, position in pipeline, corresponding m file 
   enumeration
      classification          (1)  % bigFwdPass
      myelinFix               (2)  % Myelin.runFix
      segmentation            (3)  % miniSegmentation
      correspondence          (4)  % correspondenceFinder
      overlapRemoval          (5)  % removeOverlaps
      globalSegmentID         (6)  % globalizeSegmentation
      buildSegmentMetaData    (7)  % buildSegmentMetaData
      globalCorrespondences   (8)  % globalizeCorrespondences
      graphConstruction       (9)  % graphConstruction
      SynapseDetection        (10) % SynEM.Seg.pipelineRun
      rawFeatures             (11) % connectEM.calculateRawFeatures
      classFeatures           (12) % connectEM.calculateClassFeatures
  neuriteContinuityPrediction (13) % connectEM.predictDataset
      saveGlobalSvgData       (14) % collectSvgDataOnCluster
      globalGraphStruct       (15) % collectGraphStructOnCluster
   end
end