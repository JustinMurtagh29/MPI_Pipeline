% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
wholeCellFile = fullfile(rootDir, 'aggloState', 'wholeCells_GTAxon_08_v4.mat');
outDir = '/tmpscratch/amotta/l4/2018-07-22-whole-cell-recall';

debugDir = sprintf('%s_debug-skeletons', datestr(now, 30));
debugDir = fullfile(outDir, debugDir);

nmlDirs = ...
    connectEM.WholeCell.Data.getFile({ ...
        'border-cells_axon-dendrites-split', ...
        'center-cells_axon-dendrites-split'});

segParam = struct;
segParam.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
segParam.backend = 'wkwrap';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

% Use WKW segmentation for speed
param.seg = segParam;

wholeCells = load(wholeCellFile);
wholeCells = wholeCells.wholeCells;

wholeCellAgglos = Agglo.fromSuperAgglo(wholeCells);
wholeCellLUT = Agglo.buildLUT(maxSegId, wholeCellAgglos);

%% Find NML files
nmlFiles = cellfun(@(nmlDir) ...
    dir(fullfile(nmlDir, '*.nml')), ...
    nmlDirs, 'UniformOutput', false);
nmlFiles = cellfun(@(nmlDir, nmlFiles) ...
    fullfile(nmlDir, reshape({nmlFiles.name}, [], 1)), ...
    nmlDirs, nmlFiles, 'UniformOutput', false);
nmlFiles = cat(1, nmlFiles{:});

%% Calculate path length and errors
wholeCellT = table;
wholeCellT.axonLength = nan(size(nmlFiles(:)));
wholeCellT.dendriteLength = nan(size(nmlFiles(:)));

errorData = { ...
    'nmlFile', 'cellId', 'treeName', ...
    'isAxon', 'pathLength', 'pathLengthRecalled'};
errorData = cell2struct(cell(0, numel(errorData)), errorData, 2);

for curIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curIdx};
    curNml = slurpNml(curNmlFile);
    
    curTrees = NML.buildTreeTable(curNml);
    curComments = NML.buildCommentTable(curNml);
    
    curNodes = NML.buildNodeTable(curNml);
    curNodes.coord = curNodes.coord + 1;
    
    
    % Check if trees indicate axon
    curAxonTreeName = curTrees.id(contains( ...
        curTrees.name, 'axon', 'IgnoreCase', true));
    
    % Check if comments indicate axon
    curAxonComment = curComments.node(contains( ...
        curComments.comment, 'axon', 'IgnoreCase', true));
	curAxonComment = unique(curNodes.treeId( ...
        ismember(curNodes.id, curAxonComment)));
    
    % Find smallest tree
   [~, curAxonSize] = min(cellfun( ...
        @(n) numel(n.id), curTrees.nodes));
    curAxonSize = curTrees.id(curAxonSize);
    
    curAxonId = nan;
    if ~isempty(curAxonTreeName) && ~isempty(curAxonComment)
        assert(isequal(curAxonTreeName, curAxonComment));
        curAxonId = union(curAxonTreeName, curAxonComment);
    elseif ~isempty(curAxonTreeName) || ~isempty(curAxonComment)
        curAxonId = union(curAxonTreeName, curAxonComment);
    elseif height(curTrees) > 1
        curAxonId = curAxonSize;
    end
    
    curTrees.isAxon = curTrees.id == curAxonId;
    assert(sum(curTrees.isAxon) <= 1);
    
    
    curNodes.segIds = ...
        Skeleton.getSegmentIdsOfNodes( ...
            param, curNodes.coord, 26);
    
    % Find whole cell agglomerate
    curWholeCellId = curNodes.segIds(curNodes.segIds > 0);
    curWholeCellId = nonzeros(wholeCellLUT(curWholeCellId));
    
    if isempty(curWholeCellId)
        % NOTE(amotta): This block is reached only if we've missed to
        % reconstruct a whole cell... It's existence makes me sad.
        curWholeCellId = 0;
        curWholeCellSegIds = [];
    else
        curWholeCellId = mode(curWholeCellId);
        curWholeCellSegIds = wholeCellAgglos{curWholeCellId};
    end
    
    % Determine node recall
    curNodes.isRecalled = any(ismember( ...
        curNodes.segIds, curWholeCellSegIds), 2);
    
    % Ignore nodes that are outside segmentation
    curNodes.ignore = any(curNodes.segIds < 0, 2);
    
    for curTreeIdx = 1:height(curTrees)
        curTreeName = curTrees.name{curTreeIdx};
        curTreeIsAxon = curTrees.isAxon(curTreeIdx);
        
        curTreeNodes = curTrees.nodes{curTreeIdx}.id;
       [~, curTreeNodes] = ismember(curTreeNodes, curNodes.id);
        curTreeNodes = curNodes(curTreeNodes, :);

        curTreeEdges = table;
        curTreeEdges.edge = horzcat( ...
            curTrees.edges{curTreeIdx}.source, ...
            curTrees.edges{curTreeIdx}.target);
       [~, curTreeEdges.edge] = ismember( ...
            curTreeEdges.edge, curTreeNodes.id);

        curTreeEdges.isRecalled = all(reshape( ...
            curTreeNodes.isRecalled(curTreeEdges.edge), [], 2), 2);
        curTreeEdges.ignore = any(reshape( ...
            curTreeNodes.ignore(curTreeEdges.edge), [], 2), 2);

        curTreeEdgeLengths = ...
            curTreeNodes.coord(curTreeEdges.edge(:, 1), :) ...
          - curTreeNodes.coord(curTreeEdges.edge(:, 2), :);
        curTreeEdgeLengths = curTreeEdgeLengths .* param.raw.voxelSize;
        curTreeEdgeLengths = sqrt(sum(curTreeEdgeLengths .^ 2, 2));
        
        curTreePathLength = sum(curTreeEdgeLengths);
        
        curTreeEdgeRecall = curTreeEdgeLengths;
        curTreeEdgeRecall(curTreeEdges.ignore) = 0;
        curTreeEdgeRecall = ...
            sum(curTreeEdgeRecall(curTreeEdges.isRecalled)) ...
          / sum(curTreeEdgeRecall);
        
        % Build output
        curTreeErrorData = errorData([]);
        curTreeErrorData(1).nmlFile = curNmlFile;
        curTreeErrorData(1).cellId = curWholeCellId;
        curTreeErrorData(1).treeName = curTreeName;
        curTreeErrorData(1).isAxon = curTreeIsAxon;
        curTreeErrorData(1).pathLength = curTreePathLength;
        curTreeErrorData(1).pathLengthRecalled = curTreeEdgeRecall;
        errorData = cat(1, errorData, curTreeErrorData); %#ok
    end
end

%% Total
totalAxonLengthMm = sum(wholeCellT.axonLength) / 1E6;
totalDendriteLengthMm = sum(wholeCellT.dendriteLength) / 1E6;

fprintf('Total axonal path length: %.2f mm\n', totalAxonLengthMm);
fprintf('Total dendritic path length: %.2f mm\n', totalDendriteLengthMm);
