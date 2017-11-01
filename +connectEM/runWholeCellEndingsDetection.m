function runWholeCellEndingsDetection(param,suffix)

[graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);
 
graphInput = struct();
graphInput.graph = graph;
graphInput.segmentMeta = segmentMeta;
graphInput.borderMeta = borderMeta;
graphInput.globalSegmentPCA = globalSegmentPCA;

suffix = '04'
state = '11'

connectEM.generateEndingInputDataBorderWholeCells(param,suffix,state,graphInput)
connectEM.generateEndingsBorderWholeCells(param,suffix)
connectEM.generateQueriesOfBorderWholeCells(param,suffix,state,11,graphInput)

connectEM.generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
connectEM.generateDendriteEndingsWholeCells(param,suffix)
connectEM.generateDendriteQueriesOfWholeCells(param,suffix,graphInput,8)
