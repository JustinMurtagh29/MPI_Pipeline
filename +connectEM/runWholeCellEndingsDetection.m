function runWholeCellEndingsDetection(param,suffix)

[graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);
 
graphInput = struct();
graphInput.graph = graph;
graphInput.segmentMeta = segmentMeta;
graphInput.borderMeta = borderMeta;
graphInput.globalSegmentPCA = globalSegmentPCA;
   
generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
