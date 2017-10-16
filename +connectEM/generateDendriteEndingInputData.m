function generateDendriteEndingInputData(param,stateFile,suffix)
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
    if ~exist('stateFile','var')
        stateFile = 'dendrites_03_v2_splitmerged.mat';
    end
    if ~exist('suffix','var')
        suffix = '';
    end
    % Directory with input / output data
    dataDir = fullfile(param.saveFolder, 'aggloState');
    intermediateFile = fullfile(dataDir, sprintf('dendriteDirectionality%s.mat',suffix));
    outFile = fullfile(dataDir, sprintf('dendriteEndingInputData%s.mat',suffix));

    % Load data
   [graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);
   
    % Load state of axon agglomeration and load big indices as loaded by Kevin
    dendrites = load(fullfile(dataDir, stateFile));
%     dendriteIds = find(dendrites.indBigDends);
%     dendrites = dendrites.dendrites(dendrites.indBigDends);
    dendriteIds = [1:numel(dendrites.dendrites)]';
    dendrites = dendrites.dendrites;

    % Convert to old-school agglomerates
    dendriteAgglos = arrayfun( ...
        @Agglo.fromSuperAgglo, dendrites, 'UniformOutput', false);

    % Calculate dendrite directionality
    options = struct;
    options.voxelSize = param.raw.voxelSize;
    options.bboxDist = bboxDist;
    % Just calculate for the larger 5um denrites since there occurs some issue with the smaller ones    
    directionality = connectEM.calculateDirectionalityOfAgglomerates( ...
        dendriteAgglos, graph, segmentMeta, borderMeta, globalSegmentPCA, options);
    save(intermediateFile, 'directionality', '-v7.3');

    % Save results
    out = struct;
    out.dendrites = dendrites;
    out.dendriteIds = dendriteIds;
%     out.directionality = structfun( ...
%         @(x) x(dendriteIds), directionality, 'UniformOutput', false);
    out.directionality = directionality;
    out.gitInfo = Util.gitInfo();

    Util.saveStruct(outFile, out);

end

