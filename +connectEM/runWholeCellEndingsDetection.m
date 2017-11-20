function runWholeCellEndingsDetection(param,suffix)

[graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);
 
graphInput = struct();
graphInput.graph = graph;
graphInput.segmentMeta = segmentMeta;
graphInput.borderMeta = borderMeta;
graphInput.globalSegmentPCA = globalSegmentPCA;

suffix = '05'
state = '12'

connectEM.generateEndingInputDataBorderWholeCells(param,suffix,state,graphInput)
connectEM.generateEndingsBorderWholeCells(param,suffix)
connectEM.generateQueriesOfBorderWholeCells(param,suffix,state,12,graphInput)

suffix = '04'
runID = 9

connectEM.generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
connectEM.generateDendriteEndingsWholeCells(param,suffix)
connectEM.generateDendriteQueriesOfWholeCells(param,suffix,graphInput,runID)
