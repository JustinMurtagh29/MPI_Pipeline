% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
wholeCellFile = fullfile(rootDir, 'aggloState', 'wholeCells_GTAxon_08_v4_splitHealed_v1.mat');

outputDir = '';
runId = datestr(now, 30);

searchRadius = 5E3;
boxMargin = 5E3;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

box = param.bbox .* param.raw.voxelSize(:);
box = transpose(box + [+1, -1] * boxMargin);

wholeCells = load(wholeCellFile);
wholeCells = wholeCells.wholeCells(:);

%% Detect endings
clear cur*;
for curCellId = 1:numel(wholeCells)
    curWholeCell = wholeCells(curCellId);
    curWholeCell = SuperAgglo.clean(curWholeCell);
    
    curGraph = graph(curWholeCell.edges(:, 1), curWholeCell.edges(:, 2));
    curNeighLUT = Graph.edges2Neighbors(curGraph.Edges.EndNodes);
    
    curNodes = curWholeCell.nodes(:, 1:3);
    curNodes = curNodes .* param.raw.voxelSize;
    curKd = KDTreeSearcher(curNodes);
    
    curEnding = false(size(curNodes, 1), 1);
    for curNodeId = 1:size(curNodes, 1)
        curNode = curNodes(curNodeId, :);
        if any(curNode < box(1, :)); continue; end
        if any(curNode > box(2, :)); continue; end

        % Find neighbouring node ids
        curNeighIds = rangesearch(curKd, curNode, searchRadius);
        curNeighIds = union(curNeighIds{1}, curNeighLUT{curNodeId});

        curSubGraph = curGraph.subgraph(curNeighIds);
        curSubCompIds = conncomp(curSubGraph);

        curRelNodeMask = curNeighIds == curNodeId;
        curSubCompMask = curSubCompIds == curSubCompIds(curRelNodeMask);
        curNeighIds = curNeighIds(curSubCompMask);

        % Determine if ending or not
        curPoints = curNodes(curNeighIds, :);
        curPoints = curPoints - curNode;
        curPca = pca(curPoints, 'NumComponents', 1);

        curPoints = curPoints(curNeighIds ~= curNodeId, :);
        curPoints = curPoints ./ sqrt(sum(curPoints .^ 2, 2));
        curDirs = curPoints * curPca;

        % I'd like the future to be positive
        if median(curDirs) > 0; curDirs = -curDirs; end
        if any(curDirs > 0.1); continue; end
        curEnding(curNodeId) = true;
    end
    
    if ~any(curEnding); continue; end
    
    curSkel = skeleton();
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel = curSkel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));
    
    curSkel = connectEM.Tweak.buildNml(curWholeCell, curCellId, curSkel);
    curSkel = curSkel.addBranchpoint(find(curEnding)); %#ok
    
    curSkelFileName = fullfile(outputDir, sprintf( ...
        'whole-cell-%d_run-%s.nml', curCellId, runId));
    curSkel.write(curSkelFileName);
end
