% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
wholeCellFile = fullfile(rootDir, 'aggloState', 'wholeCells_GTAxon_08_v4.mat');
dendriteFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');

probThresh = 0.75;

annotationDir = connectEM.WholeCell.Data.getFile('split-resolution');
outputDir = '';
runId = datestr(now, 30);

info = Util.runInfo();
Util.showRunInfo(info);

%% Utility function
first = @(v) v(1);
wmean = @(v, w) sum(v .* (w / sum(w)));

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

graph = Graph.load(rootDir);
graph.id = reshape(1:height(graph), [], 1);

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

curMask = graph.borderIdx > 0;
graph.borderArea = nan(height(graph), 1);
graph.borderArea(curMask) = borderAreas(graph.borderIdx(curMask));
clear borderAreas;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

dendrites = load(dendriteFile);

dendriteIds = find( ...
    dendrites.indBigDends ...
 & ~dendrites.indWholeCells ...
 & ~dendrites.indAIS);

dendrites = dendrites.dendrites;
dendAgglos = Agglo.fromSuperAgglo(dendrites);
dendLUT = Agglo.buildLUT(maxSegId, dendAgglos(dendriteIds), dendriteIds);

wholeCells = load(wholeCellFile);
wholeCells = wholeCells.wholeCells(:);

%% Generate skeletons for annotation
clear cur*;

wholeCellsPost = struct( ...
    'nodes', [], 'edges', []);

for curCellId = 1:numel(wholeCells)
    % Find and load annotations from previous rounds
    curNmlFiles = sprintf('whole-cell-%d_run-*.nml', curCellId);
    curNmlFiles = dir(fullfile(annotationDir, curNmlFiles));
    curNmlFiles = {curNmlFiles(~[curNmlFiles.isdir]).name};
    
    curAnnT = {'dendId', 'edgeId', 'tag'};
    curAnnT = table({}, {}, {}, 'VariableNames', curAnnT);
    
    for curIdx = 1:numel(curNmlFiles)
        curNmlFile = curNmlFiles{curIdx};
        curNmlFile = fullfile(annotationDir, curNmlFile);
        curNml = slurpNml(curNmlFile);
        
        % Find annotated dendrite agglomerates
        curNmlAnnT = [ ...
            'Dendrite (?<dendId>\d+)\. ', ...
            'Edge (?<edgeId>\d+) ', ...
            '.*\((?<tag>\w+)\)$'];
        curNmlAnnT = regexpi( ...
            curNml.things.name, curNmlAnnT, 'names', 'once');
        curNmlAnnT(cellfun(@isempty, curNmlAnnT)) = [];
        curNmlAnnT = struct2table(cat(1, curNmlAnnT{:}));
        
        curAnnT = cat(1, curAnnT, curNmlAnnT);
    end
    
    curAnnT.dendId = cellfun(@str2double, curAnnT.dendId);
    curAnnT.edgeId = cellfun(@str2double, curAnnT.edgeId);
    curAnnT.merge = strcmpi(curAnnT.tag, 'yes');
    
    curWholeCell = wholeCells(curCellId);
    curWholeCell = rmfield(curWholeCell, setdiff( ...
        fieldnames(curWholeCell), fieldnames(wholeCellsPost)));
    curWholeCell = SuperAgglo.clean(curWholeCell);
    
    for curIdx = 1:height(curAnnT)
        if ~curAnnT.merge(curIdx); continue; end
        
        % Patch dendrite agglomerate into whole cell
        curDend = dendrites(curAnnT.dendId(curIdx));
        curEdge = double(graph.edges(curAnnT.edgeId(curIdx), :));
        
       [curMask, curNodeId] = ismember(curEdge, curDend.nodes(:, end));
       
        curDend.nodes = cat(1, curDend.nodes, [ ...
            segPoints(curEdge(~curMask), :), curEdge(~curMask)]);
        curDend.edges = cat(1, curDend.edges, ...
            [curNodeId(curMask), size(curDend.nodes, 1)]);
        curWholeCell = SuperAgglo.merge(curWholeCell, curDend);
    end
    
    % Generate next round of annotations
    wholeCellsPost(curCellId) = curWholeCell;
    curSegIds = Agglo.fromSuperAgglo(curWholeCell);
    curLUT = logical(Agglo.buildLUT(maxSegId, {curSegIds}));

    curGraph = graph;
    curGraph.cellMask = curLUT(curGraph.edges);
    curGraph = curGraph(xor( ...
        curGraph.cellMask(:, 1), ...
        curGraph.cellMask(:, 2)), :);

    curNeighDendIds = curGraph.edges(~curGraph.cellMask);
    curNeighDendIds = setdiff(dendLUT(curNeighDendIds), 0);
    curNeighDendIds = setdiff(curNeighDendIds, curAnnT.dendId);

    curNeighT = table;
    curNeighT.dendId = curNeighDendIds(:);
    curNeighT.edgeIds = arrayfun( ...
        @(id) curGraph.id(any(ismember( ...
            curGraph.edges, dendAgglos{id}), 2)), ...
        curNeighDendIds, 'UniformOutput', false);

    curNeighT.meanEdgeProb = cellfun(@(ids) ...
        wmean(graph.prob(ids), graph.borderArea(ids)), curNeighT.edgeIds);

    curNeighT = sortrows(curNeighT, 'meanEdgeProb', 'descend');
    curNeighT(curNeighT.meanEdgeProb < probThresh, :) = [];
    if isempty(curNeighT); continue; end

    curNeighT.maxEdgeId = cellfun(@(ids) first( ...
        Util.sortBy(ids, graph.prob(ids), 'descend')), curNeighT.edgeIds);
    
    curSkel = skeleton();
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel = curSkel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));

    curSkel = Superagglos.toSkel(curWholeCell, curSkel);
    curSkel.names{end} = sprintf('Whole cell %d', curCellId);
    curSkel.colors{end} = [0.000, 0.447, 0.741, 1];

    curNumDigits = ceil(log10(1 + height(curNeighT)));
    for curDendIdx = 1:height(curNeighT)
        curDendId = curNeighT.dendId(curDendIdx);
        curEdgeId = curNeighT.maxEdgeId(curDendIdx);
        curEdge = graph.edges(curEdgeId, :);
        
        curDendAgglo = dendrites(curDendId);
       [curEdgeMask, curEdgeNodeId] = ismember( ...
            curEdge, curDendAgglo.nodes(:, end));
        
        curWcNodeId = curEdgeNodeId(curEdgeMask);
        curDendSegId = curEdge(~curEdgeMask);
        
        curDendAgglo.nodes = cat(1, curDendAgglo.nodes, ...
            [segPoints(curDendSegId, :), curDendSegId]);
        curDendAgglo.edges = cat(1, curDendAgglo.edges, ...
            [curWcNodeId, size(curDendAgglo.nodes, 1)]);
        curDendAgglo = SuperAgglo.clean(curDendAgglo, false);

        curSkel = Superagglos.toSkel(curDendAgglo, curSkel);
        curSkel = curSkel.addBranchpoint(curSkel.largestID);
        curSkel.names{end} = sprintf( ...
            '%0*d. Dendrite %d. Edge %d', ...
            curNumDigits, curDendIdx, curDendId, curEdgeId);
        curSkel.colors{end} = [0.929, 0.694, 0.1250, 1];
    end
    
    curSkelFileName = fullfile(outputDir, sprintf( ...
        'whole-cell-%d_run-%s.nml', curCellId, runId));
    curSkel.write(curSkelFileName);
end

%% Export current state of whole cells
clear cur*;
wholeCellsPost = wholeCellsPost(:);

% Mark axonal nodes as such
curAxonMasks = arrayfun(@(pre, post) ismember( ...
    post.nodes(:, end), pre.nodes(pre.axon > 0, end)), ...
    wholeCells, wholeCellsPost, 'UniformOutput', false);
[wholeCellsPost.axon] = deal(curAxonMasks{:});

out = struct;
out.info = info;
out.wholeCells = wholeCellsPost;

outFile = fullfile(outputDir, ...
    sprintf('%s_whole-cells.mat', runId));
Util.saveStruct(outFile, out);
