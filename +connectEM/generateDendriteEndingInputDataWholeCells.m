function generateDendriteEndingInputDataWholeCells(param,suffix,graphInput)
    % generateAxonEndingInputData(param)
    %   Generates a `endingInputData.mat` file in the pipeline directory
    %   which contains all the data necessary for the ending detection.
    %
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    % Parameters
    bboxDist = 1000; % in nm
    if ~exist('suffix','var')
        suffix = '';
    end
    
    if nargin < 3
        [graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
            connectEM.loadAllSegmentationData(param);
    end
       
    % Directory with input / output data
    dataDir = fullfile(param.saveFolder, 'aggloState');
    intermediateFile = fullfile(dataDir, sprintf('wholeCellsDirectionality_%s.mat',suffix));
    outFile = fullfile(dataDir, sprintf('wholeCellsEndingInputData_%s.mat',suffix));

    % Load data
    graph = graphInput.graph;
    segmentMeta = graphInput.segmentMeta;
    borderMeta = graphInput.borderMeta;
    globalSegmentPCA = graphInput.globalSegmentPCA;
   
   
    % Load state of axon agglomeration and load big indices as loaded by Kevin
    wholeCells = load(fullfile(dataDir, sprintf('wholeCells_%s.mat',suffix)));
    wholeCells = wholeCells.wholeCells;
    
    % Convert to old-school agglomerates
    wholeCellAgglos = arrayfun( ...
        @Agglo.fromSuperAgglo, wholeCells, 'UniformOutput', false);

    % Calculate dendrite directionality
    options = struct;
    options.voxelSize = param.raw.voxelSize;
    options.bboxDist = bboxDist;
    % Just calculate for the larger 5um denrites since there occurs some issue with the smaller ones    
    directionality = connectEM.calculateDirectionalityOfAgglomerates( ...
        wholeCellAgglos, graph, segmentMeta, borderMeta, globalSegmentPCA, options);
    
    % Save results
    out = struct;
    out.wholeCells = wholeCells;
    out.directionality = directionality;
    out.gitInfo = Util.gitInfo();

    Util.saveStruct(outFile, out);

end

