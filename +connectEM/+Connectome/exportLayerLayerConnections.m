% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputMapFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20190117T143833_results.mat';

dendriteFile = '/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/dendrites_wholeCells_03_v2.mat';
wholeCellFile = '/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/dendrites_wholeCells_02_v3_auto-and-manual.mat';

l4ConnRunId = '20190221T112510';
outDir = '/home/amotta/Desktop';

info = Util.runInfo();
Util.showRunInfo(info);

%% Load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);
segSizes = Seg.Global.getSegToSizeMap(param);

curData = load(outputMapFile);
axonT = struct2table(curData.axonData, 'AsArray', true);

conn = curData.info.param.connFile;
[conn, syn] = connectEM.Connectome.load(param, conn);

dendrites = load(dendriteFile);
dendrites = dendrites.dendrites(dendrites.indBigDends | dendrites.indAIS);
assert(isequal(numel(dendrites), numel(conn.dendrites)));

curData = load(wholeCellFile);
curData = curData.dendrites(curData.indWholeCells);
dendrites = cat(1, dendrites(:), curData(:));

[curDir, curFile] = fileparts(outputMapFile);
curAsiFile = sprintf('%s__%s_connectome.mat', curFile, l4ConnRunId);
curData = load(fullfile(curDir, curAsiFile));

l4SynT = curData.synT;
l4AsiT = curData.asiT;

l4AsiT = l4AsiT(l4AsiT.area > 0, :);
l4AsiT = connectEM.Consistency.Calibration.apply(l4AsiT);

%% Prepare synapse table
l4SynT.targetClass = conn.denMeta.targetClass(l4SynT.postAggloId);

l4SynT.type(:) = categorical({'Shaft'});
l4SynT.asiArea(:) = nan;

[~, curIds] = ismember(l4AsiT.id, l4SynT.id);
l4SynT.type(curIds) = l4AsiT.type;
l4SynT.asiArea(curIds) = l4AsiT.area;

%% Define utilities
wmean = @(v, w) sum(v .* (w / sum(w)), 1);
segCom = @(ids) wmean(segPoints(ids, :), segSizes(ids));

%% Process axon NML files
axonT.somaPos = nan(height(axonT), 3);
axonT.cellId = nan(height(axonT), 1);

axonT.nodes(:) = {nan(0, 3)};
axonT.edges(:) = {nan(0, 2)};

for curId = 1:height(axonT)
    curNmlFile = axonT.nmlFile{curId};
    curSkel = skeleton(curNmlFile);
    
    curAxonTreeId = curSkel.getTreeWithName('Axon');
    curDendTreeId = curSkel.getTreeWithName('Dendrite');
    
    assert(isscalar(curAxonTreeId));
    assert(isscalar(curDendTreeId));
    
    curSomaNodeId = ...
        curSkel.getNodesWithComment( ...
            'Soma', curDendTreeId, 'insensitive');
    
	axonT.somaPos(curId, :) = ...
        curSkel.nodes{curDendTreeId}(curSomaNodeId, 1:3);
    
    axonT.nodes{curId} = curSkel.nodes{curAxonTreeId}(:, 1:3);
    axonT.edges{curId} = curSkel.edges{curAxonTreeId};
end

%% Find cell IDs for axon tracings
% NOTE(amotta): We will use these coordinates in conjuction with the soma
% nodes of the axon tracings to match tracings to cell IDs. This is so
% nasty that it hurts a bit...
clear cur*;

curMask = conn.denMeta.targetClass == 'Somata';
curSomaT = conn.denMeta(curMask, {'id', 'cellId'});
curSomaT = sortrows(curSomaT, 'cellId');

curSomaT.pos = conn.dendrites(curSomaT.id);
curSomaT.pos = cellfun(segCom, curSomaT.pos, 'UniformOutput', false);
curSomaT.pos = round(cell2mat(curSomaT.pos));

[~, curIds] = pdist2( ...
    curSomaT.pos .* param.raw.voxelSize, ...
    axonT.somaPos .* param.raw.voxelSize, ...
    'squaredeuclidean', 'Smallest', 1);

axonT.cellId = curSomaT.cellId(curIds);
assert(all(axonT.cellId));

%% Load list of verified synapse IDs
clear cur*;

curAnnDir = fullfile( ...
    fileparts(fileparts(mfilename('fullpath'))), ...
    '+Consistency', '+Manual', 'annotations');

curAnnFiles = dir(fullfile(curAnnDir, 'l4-*-spine-synapses.nml'));
curAnnFiles = fullfile(curAnnDir, {curAnnFiles.name});

l4SynT.checked = false(height(l4SynT), 1);

for curId = 1:numel(curAnnFiles)
    curAnnFile = curAnnFiles{curId};
    curNml = slurpNml(curAnnFile);
    curTreeNames = curNml.things.name;
    
    curSynIds = regexp(curTreeNames, '^Syn (\d+)', 'tokens', 'once');
    curSynIds = cellfun(@str2double, cat(1, curSynIds{:}));
    
    l4SynT.checked(ismember(l4SynT.id, curSynIds)) = true;
end

%% Generate skeletons
clear cur*;

% Export synapses as
% * mst ⇒ minimum spanning tree
% * com ⇒ center of mass
curSynExportType = 'mst';
curAxonAggloExportType = '';
curDendAggloExportType = '';

assert(ismember(curSynExportType, {'', 'mst', 'com'}));
assert(ismember(curAxonAggloExportType, {'', 'mst'}));
assert(ismember(curDendAggloExportType, {'', 'agglo'}));

curMst = @(ids, skel) Skeleton.fromMST( ...
    segPoints(ids, :), param.raw.voxelSize, skel);

% Only export axons with synapses
curAxonIds = unique(l4AsiT.preAggloId);

for curAxonIdx = 1:numel(curAxonIds)
    curAxonId = curAxonIds(curAxonIdx);
    curAxonT = axonT(curAxonId, :);
    
    curAxonName = curAxonT.nmlFile{1};
   [~, curAxonName] = fileparts(curAxonName);

    curSynT = l4SynT(l4SynT.preAggloId == curAxonId, :);
    curNumDigits = ceil(log10(1 + height(curSynT)));
    if isempty(curSynT); continue; end
    
    curSkel = skeleton();
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);
    
    curSkel = curSkel.addTree( ...
        'Axon Tracing', curAxonT.nodes{1}, curAxonT.edges{1});

    % Add axon agglomerate
    if strcmpi(curAxonAggloExportType, 'mst')
        curSkel = curMst(curAxonT.segIds, curSkel);
        curSkel.names{end} = 'Axon Agglomerate';
    end
    
    % Add synapses
   [curSkel, curSynGroupId] = curSkel.addGroup('Synapses');
    
    for curSynIdx = 1:height(curSynT)        
        curSyn = curSynT(curSynIdx, :);        
        curPreSegIds = syn.synapses.presynId{curSyn.id};
        curPostSegIds = syn.synapses.postsynId{curSyn.id};
        
        if strcmpi(curSynExportType, 'mst')
            curSegIds = unique(cat(1, curPreSegIds, curPostSegIds)); 
            curSkel = curMst(curSegIds, curSkel);
        elseif strcmpi(curSynExportType, 'com')
            curPos = cat(1, segCom(curPreSegIds), segCom(curPostSegIds));
            curSkel = curSkel.addNodesAsTrees(round(mean(curPos, 1)));
        else
            continue;
        end
        
        curSkel.names{end} = sprintf( ...
            '%0*d. Synapse %d. %s onto %s %d. %f µm²', ...
            curNumDigits, curSynIdx, curSyn.id, curSyn.type, ...
            curSyn.targetClass, curSyn.postAggloId, curSyn.asiArea);
        curSkel.groupId(end) = curSynGroupId;
    end
    
    % Add presynaptic cell and postsynaptic targets
    curWcMask = ismember(curSynT.targetClass, {'Somata', 'WholeCell'});
    
    curWcIds = curSynT.checked & curWcMask;
    curWcIds = conn.denMeta.cellId(curSynT.postAggloId(curWcIds));
    curWcIds = union(curAxonT.cellId, curWcIds, 'stable');
    
    curWcNames = arrayfun( ...
        @(id) sprintf('Cell %d', id), ...
        curWcIds, 'UniformOutput', false);
    curWcNames{1} = strcat(curWcNames{1}, ' (with Axon)');
    
    curDendIds = curSynT.checked & not(curWcMask);
    curDendIds = unique(conn.denMeta.id(curSynT.postAggloId(curDendIds)));
    
    curDendNames = arrayfun( ...
        @(target, id) sprintf('%s %d', target, id), ...
        conn.denMeta.targetClass(curDendIds), curDendIds, ...
        'UniformOutput', false);
    
    % NOTE(amotta): The whole cells were appended to the dendrites.
    curDendIds = cat(1, curWcIds(:) + numel(conn.dendrites), curDendIds(:));
    curDendNames = cat(1, curWcNames(:), curDendNames(:));
    
    for curDendIdx = 1:numel(curDendIds)
        if isempty(curDendAggloExportType); continue; end
        
        curDendName = curDendNames(curDendIdx);
        curDendId = curDendIds(curDendIdx);
        
        curDend = dendrites(curDendId);
        curDend = SuperAgglo.clean(curDend);
        
       [curAgglos, curSkels] = ...
            Superagglos.splitIntoAgglosAndFlights(curDend);
       [curSkel, curGroupId] = curSkel.addGroup(curDendName);
        
       [curSkel, curAggloGroupId] = ...
            curSkel.addGroup('Agglomerates', curGroupId);
        for curWcAgglo = reshape(curAgglos{1}, 1, [])
            curSkel = curSkel.addTree( ...
                'Agglomerate', curWcAgglo.nodes, curWcAgglo.edges);
            curSkel.groupId(end) = curAggloGroupId;
        end
        
       [curSkel, curSkelGroupId] = ...
            curSkel.addGroup('Tracings', curGroupId);
        for curWcSkel = reshape(curSkels{1}, 1, [])
            curSkel = curSkel.addTree( ...
                'Tracing', curWcSkel.nodes, curWcSkel.edges);
            curSkel.groupId(end) = curSkelGroupId;
        end
    end
    
    curOutFile = fullfile(outDir, strcat(curAxonName, '.nml'));
    curSkel.write(curOutFile);
end
