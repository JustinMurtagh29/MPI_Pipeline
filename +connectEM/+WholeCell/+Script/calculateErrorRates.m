% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

nmlDirs = ...
    connectEM.WholeCell.Data.getFile({ ...
        'border-cells_axon-dendrites-split', ...
        'center-cells_axon-dendrites-split'});

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% Find NML files
nmlFiles = cellfun(@(nmlDir) ...
    dir(fullfile(nmlDir, '*.nml')), ...
    nmlDirs, 'UniformOutput', false);
nmlFiles = cellfun(@(nmlDir, nmlFiles) ...
    fullfile(nmlDir, reshape({nmlFiles.name}, [], 1)), ...
    nmlDirs, nmlFiles, 'UniformOutput', false);
nmlFiles = cat(1, nmlFiles{:});

%% Calculate path length
wholeCellT = table;
wholeCellT.axonLength = nan(size(nmlFiles(:)));
wholeCellT.dendriteLength = nan(size(nmlFiles(:)));

for curIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curIdx};
    curNml = slurpNml(curNmlFile);
    
    curNodes = NML.buildNodeTable(curNml);
    curNodes.coord = curNodes.coord .* param.raw.voxelSize;
    
    curTrees = NML.buildTreeTable(curNml);
    curComments = NML.buildCommentTable(curNml);
    
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
    
    curAxonId = [];
    if ~isempty(curAxonTreeName) && ~isempty(curAxonComment)
        assert(isequal(curAxonTreeName, curAxonComment));
        curAxonId = union(curAxonTreeName, curAxonComment);
    elseif ~isempty(curAxonTreeName) || ~isempty(curAxonComment)
        curAxonId = union(curAxonTreeName, curAxonComment);
    elseif height(curTrees) > 1
        curAxonId = curAxonSize;
    end
    
   	curAxonId = find(curTrees.id == curAxonId);
    assert(numel(curAxonId) <= 1);
    
    curTrees.isAxon = false(height(curTrees), 1);
    curTrees.isAxon(curAxonId) = true;
    
    % Calculate path length
    curTrees.edges = cellfun( ...
        @(e) horzcat(e.source, e.target), ...
        curTrees.edges, 'UniformOutput', false);
    
    curTrees.pathLength = nan(height(curTrees), 1);
    for curTreeIdx = 1:height(curTrees)
        curEdges = curTrees.edges{curTreeIdx};
       [~, curEdges] = ismember(curEdges, curNodes.id);
        
        curPathLength = ...
            curNodes.coord(curEdges(:, 1), :) ...
          - curNodes.coord(curEdges(:, 2), :);
        curPathLength = sum(sqrt(sum(curPathLength .^ 2, 2)));
        curTrees.pathLength(curTreeIdx) = curPathLength;
    end
    
    curTotalLength = sum(curTrees.pathLength);
    curAxonLength = sum(curTrees.pathLength(curAxonId));
    wholeCellT.axonLength(curIdx) = curAxonLength;
    wholeCellT.dendriteLength(curIdx) = curTotalLength - curAxonLength;
end

%% Total
totalAxonLengthMm = sum(wholeCellT.axonLength) / 1E6;
totalDendriteLengthMm = sum(wholeCellT.dendriteLength) / 1E6;

fprintf('Total axonal path length: %.2f mm\n', totalAxonLengthMm);
fprintf('Total dendritic path length: %.2f mm\n', totalDendriteLengthMm);
