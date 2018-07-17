% ================================================================================
% PipelineStep class - allows for setting a specific starting point in the pipeline
% if for example some steps were already done before like in the case of an error mid-pipeline 
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
% ================================================================================
classdef PipelineStep < uint8

   % <step name corresponding to job name> (position in pipeline) % corresponding m file                                                                                     
   enumeration
      Classification          (1)  % bigFwdPass
      MyelinDetection         (2)  % bigFwdPassMyelinCodat
      MyelinFix               (3)  % Myelin.runFix
      Segmentation            (4)  % miniSegmentation
      OverlapRemoval          (5)  % removeOverlaps
      GlobalSegmentID         (6)  % globalizeSegmentation
      Correspondence          (7)  % correspondenceFinderGlobal
      BuildSegmentMetaData    (8)  % buildSegmentMetaData
      SegmentationPyramid     (9)  % createResolutionPyramid
      CompressSegmentation    (10)  % compressSegmentation
      GraphConstruction       (11) % graphConstruction
      SynapseDetection        (12) % SynEM.Seg.pipelineRun
      RawFeatures             (13) % connectEM.calculateRawFeatures
      ClassFeatures           (14) % connectEM.calculateClassFeatures
  NeuriteContinuityPrediction (15) % connectEM.predictDataset
      SaveGlobalSvgData       (16) % collectSvgDataOnCluster
      GlobalGraphStruct       (17) % collectGraphStructOnCluster
      
      % additional steps from comment at bottom of pipeline
      HeuristicLookup         (18) % connectEM.getHeuristicResult
      Agglomeration           (19) % connectEM.agglomerate
   end
end

