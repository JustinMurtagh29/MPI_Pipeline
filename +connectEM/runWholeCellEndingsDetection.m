function runWholeCellEndingsDetection(param,suffix)

[graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);
 
graphInput = struct();
graphInput.graph = graph;
graphInput.segmentMeta = segmentMeta;
graphInput.borderMeta = borderMeta;
graphInput.globalSegmentPCA = globalSegmentPCA;

suffix = '03'
state = '10'

connectEM.generateEndingInputDataBorderWholeCells(param,suffix,state,graphInput)
connectEM.generateEndingsBorderWholeCells(param,suffix)
connectEM.generateQueriesOfBorderWholeCells(param,suffix,state,10,graphInput)

connectEM.generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
connectEM.generateDendriteEndingsWholeCells(param,suffix)
connectEM.generateDendriteQueriesOfWholeCells(param,suffix,graphInput,8)
