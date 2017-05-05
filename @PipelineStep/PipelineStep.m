% ================================================================================
% PipelineStep class - allows for setting a specific starting point in the pipeline
% if for example some steps were already done before like in the case of an error mid-pipeline 
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
% ================================================================================
classdef PipelineStep < uint8
   
   % <step name corresponding to job name> (position in pipeline) % corresponding m file 
   enumeration
      Classification          (1)  % bigFwdPass
      MyelinFix               (2)  % Myelin.runFix
      Segmentation            (3)  % miniSegmentation
      Correspondence          (4)  % correspondenceFinder
      OverlapRemoval          (5)  % removeOverlaps
      GlobalSegmentID         (6)  % globalizeSegmentation
      BuildSegmentMetaData    (7)  % buildSegmentMetaData
      GlobalCorrespondences   (8)  % globalizeCorrespondences
      GraphConstruction       (9)  % graphConstruction
      SynapseDetection        (10) % SynEM.Seg.pipelineRun
      RawFeatures             (11) % connectEM.calculateRawFeatures
      SegmentationPyramid     (12) % createResolutionPyramid
      ClassFeatures           (13) % connectEM.calculateClassFeatures
  NeuriteContinuityPrediction (14) % connectEM.predictDataset
      SaveGlobalSvgData       (15) % collectSvgDataOnCluster
      GlobalGraphStruct       (16) % collectGraphStructOnCluster
      
      % additional steps from comment at bottom of pipeline
      HeuristicLookup         (17) % connectEM.getHeuristicResult
      Agglomeration           (18) % connectEM.agglomerate
   end
end