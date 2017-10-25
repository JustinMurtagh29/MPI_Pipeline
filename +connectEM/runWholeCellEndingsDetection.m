function runWholeCellEndingsDetection(param,suffix)

[graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);
 
graphInput = struct();
graphInput.graph = graph;
graphInput.segmentMeta = segmentMeta;
graphInput.borderMeta = borderMeta;
graphInput.globalSegmentPCA = globalSegmentPCA;

suffix = '02'

connectEM.generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
connectEM.generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
connectEM.generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
