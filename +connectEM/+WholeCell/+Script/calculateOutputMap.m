% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

runId = datestr(now, 30);
outDir = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps';

nmlDirs = ...
    connectEM.WholeCell.Data.getFile({ ...
        'border-cells_axon-dendrites-split', ...
        'center-cells_axon-dendrites-split'});

segParam = struct;
segParam.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
segParam.backend = 'wkwrap';

info = Util.runInfo();
Util.showRunInfo(info);

%% Define helper functions
wmean = @(v, w) sum(v .* (w / sum(w)));

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.seg = segParam;

[conn, syn] = connectEM.Connectome.load(param, connFile);

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

%% Preprocess synapses
dendLUT = repelem( ...
    1:numel(conn.dendrites), ...
    cellfun(@numel, conn.dendrites));
dendLUT = accumarray( ...
    cell2mat(conn.dendrites), dendLUT(:), ...
   [maxSegId, 1], @(ids) {ids(:)}, {zeros(0, 1)});

syn.synapses.dendIds = cellfun( ...
    @(ids) unique(cell2mat(dendLUT(ids))), ...
    syn.synapses.postsynId, 'UniformOutput', false);

synLUT = repelem( ...
    1:height(syn.synapses), ...
    cellfun(@numel, syn.synapses.presynId));
synLUT = accumarray( ...
    cell2mat(syn.synapses.presynId), synLUT(:), ...
   [maxSegId, 1], @(ids) {ids(:)}, {zeros(0, 1)});

%% Find NML files
nmlFiles = cellfun(@(nmlDir) ...
    dir(fullfile(nmlDir, '*.nml')), ...
    nmlDirs, 'UniformOutput', false);
nmlFiles = cellfun(@(nmlDir, nmlFiles) ...
    fullfile(nmlDir, reshape({nmlFiles.name}, [], 1)), ...
    nmlDirs, nmlFiles, 'UniformOutput', false);
nmlFiles = cat(1, nmlFiles{:});

%% For development purposes
% nmlFiles = connectEM.WholeCell.Data.getFile(fullfile( ...
% 	  'border-cells_axon-dendrites-split', '5a796fcd6700004b172d94e9.nml'));
% nmlFiles = {nmlFiles};

%% Collect output synapses
clear cur*;

preSynAgglos = syn.synapses.presynId;
axonData = cell(size(nmlFiles));

parfor curIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curIdx};
    curNml = slurpNml(curNmlFile);
    
    curTree = NML.buildTreeTable(curNml);
    
    % Validate tree names and find axon
    assert(all(ismember(curTree.name, {'Axon', 'Dendrite'})));
    curTree = curTree(strcmpi(curTree.name, 'Axon'), :);
    
    if isempty(curTree); continue; end
    assert(height(curTree) == 1);
    
    curNodes = NML.buildNodeTable(curNml);
    curNodes(curNodes.treeId ~= curTree.id, :) = [];
    curNodes.coord = curNodes.coord + 1;
    
    % Find soma nodes
    curSomaNodeIds = NML.buildCommentTable(curNml);
    curSomaNodeIds = curSomaNodeIds.node( ...
        strcmpi(curSomaNodeIds.comment, 'soma'));
    
    curNodes.isSoma = ismember(curNodes.id, curSomaNodeIds);
    curSomaNodeId = find(curNodes.isSoma);
    assert(isscalar(curSomaNodeId));
    
    curEdges = table;
    curEdges.edge = horzcat( ...
        curTree.edges{1}.source, ...
        curTree.edges{1}.target);
   [~, curEdges.edge] = ismember( ...
        curEdges.edge, curNodes.id);
    
    curEdges.length = ...
        curNodes.coord(curEdges.edge(:, 1), :) ...
      - curNodes.coord(curEdges.edge(:, 2), :);
    curEdges.length = curEdges.length .* param.raw.voxelSize; %#ok
    curEdges.length = sqrt(sum(curEdges.length .^ 2, 2));
    curPathLength = sum(curEdges.length);
    
    curGraph = graph( ...
        curEdges.edge(:, 1), curEdges.edge(:, 2), ...
        curEdges.length, height(curNodes));
    curNodes.somaDist(:) = curGraph.distances(curSomaNodeId);
    
    curNodes.segIds = ...
        Skeleton.getSegmentIdsOfNodes(param, curNodes.coord, 26);
    curAxonSegIds = unique(curNodes.segIds(curNodes.segIds > 0));
    
    curSynT = table;
    curSynT.id = cell2mat(synLUT(curAxonSegIds)); %#ok
    curSynT.somaDist = cellfun( ...
        @(segIds) wmean( ...
            curNodes.somaDist, ...
            sum(ismember(curNodes.segIds, segIds), 2)), ...
        preSynAgglos(curSynT.id)); %#ok
    curSynT = sortrows(curSynT, 'somaDist');
    
    curAxonData = struct;
    curAxonData.nmlFile = curNmlFile;
    curAxonData.pathLength = curPathLength;
    curAxonData.segIds = curAxonSegIds(:);
    curAxonData.synapses = curSynT;
    axonData{curIdx} = curAxonData;
end

axonData = cat(1, axonData{:});

%% Save results
out = struct;
out.info = info;
out.axonData = axonData;

outFile = fullfile(outDir, sprintf('%s_results.mat', runId));
Util.saveStruct(outFile, out);
Util.protect(outFile);
