function generateDendriteEndingInputDataWholeCells(param,suffix)
    % generateAxonEndingInputData(param)
    %   Generates a `endingInputData.mat` file in the pipeline directory
    %   which contains all the data necessary for the ending detection.
    %
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    % Parameters
    minDendriteLength = 5000; % in nm
    bboxDist = 1000; % in nm
    if ~exist('suffix','var')
        suffix = '';
    end
    % Directory with input / output data
    dataDir = fullfile(param.saveFolder, 'aggloState');
    intermediateFile = fullfile(dataDir, sprintf('wholeCellsDirectionality_%s.mat',suffix));
    outFile = fullfile(dataDir, sprintf('wholeCellsEndingInputData_%s.mat',suffix));

    % Load data
   [graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);
   
    % Load state of axon agglomeration and load big indices as loaded by Kevin
    wholeCells = load(fullfile(dataDir, sprintf('wholeCells_%s.mat',suffix)));
%     dendriteIds = find(dendrites.WholeCellId);
%     dendrites = dendrites.dendrites(dendrites.WholeCellId);
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
    out.dendrites = dendrites;
%     out.dendriteIds = dendriteIds;
    out.directionality = directionality;
    out.gitInfo = Util.gitInfo();

    Util.saveStruct(outFile, out);

end

