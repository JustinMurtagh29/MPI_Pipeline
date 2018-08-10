% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
wholeCellFile = fullfile(rootDir, 'aggloState', 'wholeCells_GTAxon_08_v4_splitHealed_v1.mat');

outDir = '/tmpscratch/amotta/l4/2018-07-28-whole-cell-recall';
runId = '20180726T150513';

debugDir = sprintf('%s_debug-skeletons', runId);
debugDir = fullfile(outDir, debugDir);

nmlDirs = ...
    connectEM.WholeCell.Data.getFile({ ...
        'border-cells_axon-dendrites-split', ...
        'center-cells_axon-dendrites-split'});

splitMinLenNm = 5000;
splitBoxMarginNm = 5000;

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

% Only count splits if they enter this bounding box
splitBox = round(splitBoxMarginNm ./ param.raw.voxelSize);
splitBox = param.bbox + [+1, -1] .* splitBox(:);

%% Find NML files
nmlFiles = cellfun(@(nmlDir) ...
    dir(fullfile(nmlDir, '*.nml')), ...
    nmlDirs, 'UniformOutput', false);
nmlFiles = cellfun(@(nmlDir, nmlFiles) ...
    fullfile(nmlDir, reshape({nmlFiles.name}, [], 1)), ...
    nmlDirs, nmlFiles, 'UniformOutput', false);
nmlFiles = cat(1, nmlFiles{:});

%% Prepare for debugging
if ~exist(debugDir, 'dir')
    mkdir(debugDir);
else
    % Write skeletons only once
    debugDir = [];
end

%% Calculate path length and errors
errorData = cell(numel(nmlFiles), 1);
parfor curIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curIdx};
    
    try
    curNml = slurpNml(curNmlFile);
    
    curTrees = NML.buildTreeTable(curNml);
    curComments = NML.buildCommentTable(curNml);
    
    curNodes = NML.buildNodeTable(curNml);
    curNodes.coord = curNodes.coord + 1;
    
    % Validate tree names and find axon
    assert(all(ismember(curTrees.name, {'Axon', 'Dendrite'})));
    curTrees.isAxon = strcmpi(curTrees.name, 'Axon');
    
    % Find soma nodes
    curNodes.isSoma(:) = false;
   [~, curRowIds] = ismember(curComments.node, curNodes.id);
    curNodes.isSoma(curRowIds) = strcmpi(curComments.comment, 'soma');
    
    curNodes.segIds = ...
        Skeleton.getSegmentIdsOfNodes( ...
            param, curNodes.coord, 26);
    
    % Find whole cell agglomerate
    curWholeCellId = curNodes.segIds(curNodes.segIds > 0);
    curWholeCellId = nonzeros(wholeCellLUT(curWholeCellId)); %#ok
    
    if isempty(curWholeCellId)
        % NOTE(amotta): This block is reached only if we've missed to
        % reconstruct a whole cell... It's existence makes me sad.
        curWholeCellId = 0;
        curWholeCellSegIds = [];
    else
        curWholeCellId = mode(curWholeCellId);
        curWholeCellSegIds = wholeCellAgglos{curWholeCellId}; %#ok
    end
    
    % Determine node recall
    curNodes.isRecalled = any(ismember( ...
        curNodes.segIds, curWholeCellSegIds), 2);
    
    % Ignore nodes that are outside segmentation
    curNodes.ignore = any(curNodes.segIds < 0, 2);
    
    
    % Prepare debug skeleton
    curDebugSkel = [];
    if ~isempty(debugDir) && curWholeCellId ~= 0
        curDebugSkel = skeleton();
        curDebugSkel = Skeleton.setParams4Pipeline(curDebugSkel, param);
        curDebugSkel = curDebugSkel.setDescription(sprintf( ...
            '%s (%s)', info.filename, info.git_repos{1}.hash)); %#ok
        
        curDebugSkel = Superagglos.toSkel( ...
            SuperAgglo.clean(wholeCells(curWholeCellId)), curDebugSkel); %#ok
        curDebugSkel.names{end} = sprintf('Whole cell %d', curWholeCellId);
    end
    
    
    curErrorData = cell(height(curTrees), 1);
    for curTreeIdx = 1:height(curTrees)
        curTreeName = curTrees.name{curTreeIdx};
        curTreeIsAxon = curTrees.isAxon(curTreeIdx);
        
        curTreeNodes = curTrees.nodes{curTreeIdx}.id;
       [~, curTreeNodes] = ismember(curTreeNodes, curNodes.id);
        curTreeNodes = curNodes(curTreeNodes, :);
        
        curTreeSomaNodeId = find(curTreeNodes.isSoma);
        assert(isscalar(curTreeSomaNodeId));

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

        curTreeEdges.length = ...
            curTreeNodes.coord(curTreeEdges.edge(:, 1), :) ...
          - curTreeNodes.coord(curTreeEdges.edge(:, 2), :);
        curTreeEdges.length = curTreeEdges.length .* param.raw.voxelSize;
        curTreeEdges.length = sqrt(sum(curTreeEdges.length .^ 2, 2));
        
        
        % NOTE(amotta): To determine the number of splits, we count the
        % number of non-recalled connected components above a certain size
        % threshold. The size threshold aims to discard FP split detections
        % in, e.g., sporadically missed segments in mitochondria.
        %   Since nuclei are not part of the whole cells, the
        % "soma" node placed in the center of the soma / nucleus is treated
        % separately.
        %   Furthermore, we only count missed segments that reach at least
        % 5 µm into the dataset as splits. This is the criterion used by
        % MB, when he detected open endings.
        curGraph = ~( ...
            curTreeEdges.isRecalled ...
          | curTreeEdges.ignore);
        curGraph = graph( ...
            curTreeEdges.edge(curGraph, 1), ...
            curTreeEdges.edge(curGraph, 2), ...
            curTreeEdges.length(curGraph), ...
            height(curTreeNodes));
        
        curGraphComps = curGraph.conncomp();
        curGraphCompIds = unique(curGraph.Edges.EndNodes(:));
        curGraphCompIds = unique(curGraphComps(curGraphCompIds));
        
        curGraphCompLengths = @(id) sum( ...
            curGraph.subgraph(curGraphComps == id).Edges.Weight);
        curGraphCompNodes = @(var, id) ...
            curTreeNodes.(var)(curGraphComps == id, :);
        curGraphCompIsSoma = @(id) any(curGraphCompNodes('isSoma', id));
        curGraphCompIsInBox = @(id) any( ...
            all(curGraphCompNodes('coord', id) >= splitBox(:, 1)', 2) ...
          & all(curGraphCompNodes('coord', id) <= splitBox(:, 2)', 2)); %#ok
        
        curGraphCompIsSplit = arrayfun(@(id) ...
            curGraphCompLengths(id) >= splitMinLenNm ...
          & not(curGraphCompIsSoma(id)) ...
          & curGraphCompIsInBox(id), ...
            curGraphCompIds);
        
        % Find nodes corresponding to splits
        curNodeSomaDists = graph( ...
            curTreeEdges.edge(:, 1), curTreeEdges.edge(:, 2), ...
            curTreeEdges.length, height(curTreeNodes));
        curNodeSomaDists = curNodeSomaDists.distances(curTreeSomaNodeId);
        
        curFirst = @(ids) ids(1);
        curSplitNodes = @(ids) curFirst( ...
            Util.sortBy(ids, curNodeSomaDists(ids), 'ascend'));
        curSplitNodes = arrayfun( ...
            @(id) curSplitNodes(find(curGraphComps == id)), ...
            curGraphCompIds(curGraphCompIsSplit)); %#ok
        
        curSplitCount = sum(curGraphCompIsSplit);
        
        
        curTreePathLength = sum( ...
            curTreeEdges.length);
        curTreePathLengthValid = sum( ...
            curTreeEdges.length( ...
                ~curTreeEdges.ignore));
        curTreePathLengthRecalled = sum( ...
            curTreeEdges.length( ...
                ~curTreeEdges.ignore ...
               & curTreeEdges.isRecalled));
        
        
        % Build output
        curTreeErrorData = struct;
        curTreeErrorData.nmlFile = curNmlFile;
        curTreeErrorData.cellId = curWholeCellId;
        curTreeErrorData.treeName = curTreeName;
        curTreeErrorData.isAxon = curTreeIsAxon;
        curTreeErrorData.pathLength = curTreePathLength;
        curTreeErrorData.pathLengthValid = curTreePathLengthValid;
        curTreeErrorData.pathLengthRecalled = curTreePathLengthRecalled;
        curTreeErrorData.splitCount = curSplitCount;
        curErrorData{curTreeIdx} = curTreeErrorData;
        
        
        % Generate debug NML file
        if ~isempty(curDebugSkel)
            curDebugSkelBranchPointIds = curDebugSkel.largestID + ...
                find(~curTreeNodes.ignore & ~curTreeNodes.isRecalled);
            
            curDebugSkelComments = repelem({''}, height(curTreeNodes), 1);
            curDebugSkelComments(curSplitNodes) = {'Split'};
            
            curNeuriteType = {'Dendrite', 'Axon'};
            curNeuriteType = curNeuriteType{1 + curTreeIsAxon};
            curDebugSkelTreeName = sprintf( ...
                'Tracing %0*d. %s (%s)', ...
                ceil(log10(1 + height(curTrees))), ...
                curTreeIdx, curTreeName, curNeuriteType);
            
            curDebugSkel = curDebugSkel.addTree( ...
                curDebugSkelTreeName, curTreeNodes.coord, ...
                curTreeEdges.edge, [], [], curDebugSkelComments);
            
            % Mark "missed" nodes as branch points
            curDebugSkel = ...
                curDebugSkel.addBranchpoint( ...
                    curDebugSkelBranchPointIds);
        end
    end
    
    if ~isempty(curDebugSkel)
       [~, curDebugSkelFileName] = fileparts(curNmlFile);
        curNumDigits = ceil(log10(1 + numel(nmlFiles))); %#ok
        
        curDebugSkelFileName = fullfile(debugDir, sprintf( ...
            '%0*d_%s.nml', curNumDigits, curIdx, curDebugSkelFileName));
        curDebugSkel.write(curDebugSkelFileName);
    end
    
    curErrorData = vertcat(curErrorData{:});
    errorData{curIdx} = curErrorData;
    
    catch
        error('Error in %s', curNmlFile);
    end
end

errorData = vertcat(errorData{:});
errorData = struct2table(errorData);

out = struct;
out.info = info;
out.errorData = errorData;

outFile = sprintf('%s_error-data.mat', runId);
outFile = fullfile(outDir, outFile);

Util.saveStruct(outFile, out);
Util.protect(outFile);

%% Evaluation
numTracings = numel(unique(errorData.nmlFile)) %#ok
tracingsWithoutWholeCell = unique( ...
    errorData.nmlFile(~errorData.cellId)) %#ok
wholeCellsWithoutTracing = setdiff( ...
    1:numel(wholeCells), errorData.cellId) %#ok

wholeCellsWithMultipleTracings = unique(errorData( ...
    errorData.cellId > 0, {'nmlFile', 'cellId'}), 'rows');
[wholeCellsWithMultipleTracings, ~, uniIds] = ...
    unique(wholeCellsWithMultipleTracings.cellId);
wholeCellsWithMultipleTracings(accumarray(uniIds, 1) < 2) = [] %#ok

totalAxonLengthMm = sum( ...
    errorData.pathLength(errorData.isAxon)) / 1E6 %#ok
totalAxonLengthValidMm = sum( ...
    errorData.pathLengthValid(errorData.isAxon)) / 1E6 %#ok
totalAxonLengthRecalledMm = sum( ...
    errorData.pathLengthRecalled(errorData.isAxon)) / 1E6 %#ok
totalAxonRecall = ...
    totalAxonLengthRecalledMm / totalAxonLengthValidMm %#ok

totalDendriteLengthMm = sum( ...
    errorData.pathLength(~errorData.isAxon)) / 1E6 %#ok
totalDendriteLengthValidMm = sum( ...
    errorData.pathLengthValid(~errorData.isAxon)) / 1E6 %#ok
totalDendriteLengthRecalledMm = sum( ...
    errorData.pathLengthRecalled(~errorData.isAxon)) / 1E6 %#ok
totalDendriteRecall = ...
    totalDendriteLengthRecalledMm / totalDendriteLengthValidMm %#ok

%% Debug recall
clear cur*;
curIsAxon = false;

curData = errorData;
curData(curData.isAxon ~= curIsAxon, :) = [];

curData.missedPathLength = ...
    curData.pathLengthValid - curData.pathLengthRecalled;
curData = sortrows(curData, 'missedPathLength', 'descend');

fprintf('\n');
disp(curData);
