function runWholeCellEndingsDetection(param,suffix)

[graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);
 
graphInput = struct();
graphInput.graph = graph;
graphInput.segmentMeta = segmentMeta;
graphInput.borderMeta = borderMeta;
graphInput.globalSegmentPCA = globalSegmentPCA;

suffix = '03'

connectEM.generateEndingInputDataBorderWholeCells(param,suffix,graphInput)
connectEM.generateEndingsBorderWholeCells(param,suffix)
connectEM.generateQueriesOfBorderWholeCells(param,suffix,graphInput,1)

connectEM.generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
connectEM.generateDendriteEndingsWholeCells(param,suffix)
connectEM.generateDendriteQueriesOfWholeCells(param,suffix,graphInput,8)
