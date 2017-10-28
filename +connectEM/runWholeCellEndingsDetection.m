function runWholeCellEndingsDetection(param,suffix)

[graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);
 
graphInput = struct();
graphInput.graph = graph;
graphInput.segmentMeta = segmentMeta;
graphInput.borderMeta = borderMeta;
graphInput.globalSegmentPCA = globalSegmentPCA;

suffix = '02'
state = '09'

connectEM.generateEndingInputDataBorderWholeCells(param,suffix,state,graphInput)
connectEM.generateEndingsBorderWholeCells(param,suffix)
connectEM.generateQueriesOfBorderWholeCells(param,suffix,state,9,graphInput)

connectEM.generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
connectEM.generateDendriteEndingsWholeCells(param,suffix)
connectEM.generateDendriteQueriesOfWholeCells(param,suffix,graphInput,8)
