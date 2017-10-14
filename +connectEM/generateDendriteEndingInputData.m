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
        stateFile = 'dendrites_03_v2.mat';
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
%     graph=load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNewNew.mat')
% borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
% segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'point', 'maxSegId', 'cubeIdx', 'centroid', 'box');
% load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat', 'p');
% segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
% globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
% [graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
% graph.neighProb = cellfun(@(x)graph.prob(x), neighboursIdx, 'uni', 0);
% graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);

   
   
    % Load state of axon agglomeration and load big indices as loaded by Kevin
    dendrites = load(fullfile(dataDir, stateFile));
    dendriteIds = find(dendrites.indBigDends);
    dendrites = dendrites.dendrites;

    % Convert to old-school agglomerates
%     dendriteAgglos = arrayfun( ...
%         @Agglo.fromSuperAgglo, dendrites, 'UniformOutput', false);
    dendriteAgglos = arrayfun(@(x)x.nodes(:,4),dendrites,'uni',0);

    % Calculate axon directionality
    options = struct;
    options.voxelSize = param.raw.voxelSize;
    options.bboxDist = bboxDist;

    directionality = connectEM.calculateDirectionalityOfAgglomerates( ...
        dendriteAgglos, graph, segmentMeta, borderMeta, globalSegmentPCA, options);
    save(intermediateFile, 'directionality', '-v7.3');

    % Save results
    out = struct;
    out.dendrites = dendrites(dendriteIds);
    out.dendriteIds = dendriteIds;
    out.directionality = structfun( ...
        @(x) x(dendriteIds), directionality, 'UniformOutput', false);
    out.gitInfo = Util.gitInfo();

    Util.saveStruct(outFile, out);

end

