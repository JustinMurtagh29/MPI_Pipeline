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

%% Prepare for debugging
if ~isempty(debugDir)
    mkdir(debugDir);
end

%% Calculate path length and errors
errorData = cell(numel(nmlFiles), 1);
parfor curIdx = 1:numel(nmlFiles)
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
        curErrorData{curTreeIdx} = curTreeErrorData;
        
        
        % Generate debug NML file
        if ~isempty(curDebugSkel)
            curDebugSkelBranchPointIds = curDebugSkel.largestID + ...
                find(~curTreeNodes.ignore & ~curTreeNodes.isRecalled);
            
            curNeuriteType = {'Dendrite', 'Axon'};
            curNeuriteType = curNeuriteType{1 + curTreeIsAxon};
            curDebugSkelTreeName = sprintf( ...
                'Tracing %0*d. %s (%s)', ...
                ceil(log10(1 + height(curTrees))), ...
                curTreeIdx, curTreeName, curNeuriteType);
            
            curDebugSkel = curDebugSkel.addTree( ...
                curDebugSkelTreeName, ...
                curTreeNodes.coord, ...
                curTreeEdges.edge);
            
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
end

errorData = vertcat(errorData{:});
errorData = struct2table(errorData);

%% Evaluation
numTracings = numel(nmlFiles) %#ok
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
