% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191203T021242-results_20191203T021242-results-auto-spines-v2_SynapseAgglomerates--20191203T021242-results--20191203T021242-results-auto-spines-v2--v1.mat');

outVersion = 1;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = load(connFile);
syn = load(conn.info.param.synFile);
axons = Util.load(conn.info.param.axonFile, 'axons');
dendrites = Util.load(conn.info.param.dendriteFile, 'dendrites');

segSizes = Seg.Global.getSegToSizeMap(param);

%% Calculate axonal intersynapse distances
clear cur*;

curConn = conn.connectome;
curAxonIds = cellfun(@numel, curConn.synIdx);
curAxonIds = repelem(curConn.edges(:, 1), curAxonIds);
curSynIds = cell2mat(curConn.synIdx);

curSynIds = accumarray( ...
    curAxonIds, curSynIds, size(conn.axons), ...
    @(ids) {ids(:)}, {zeros(0, 1, 'like', curSynIds)});

[axonSynToSynDists, axonSynIds] = ...
    Synapse.calculateIntersynapseDistances( ...
        axons, curSynIds, syn.synapses.presynId, ...
        'voxelSize', param.raw.voxelSize, 'segWeights', segSizes);

%% Calculate dendritic intersynapse distances
clear cur*;

curConn = conn.connectome;
curDendIds = cellfun(@numel, curConn.synIdx);
curDendIds = repelem(curConn.edges(:, 2), curDendIds);
curSynIds = cell2mat(curConn.synIdx);

curSynIds = accumarray( ...
    curDendIds, curSynIds, size(conn.dendrites), ...
    @(ids) {ids(:)}, {zeros(0, 1, 'like', curSynIds)});

[dendSynToSynDists, dendSynIds] = ...
    Synapse.calculateIntersynapseDistances( ...
        dendrites, curSynIds, syn.synapses.postsynId, ...
        'voxelSize', param.raw.voxelSize, 'segWeights', segSizes);
    
%% Save result
clear cur*;

out = struct;
out.info = info;
out.axonSynToSynDists = axonSynToSynDists;
out.axonSynIds = axonSynIds;
out.dendSynToSynDists = dendSynToSynDists;
out.dendSynIds = dendSynIds;

outFile = sprintf('__intersynapse_v%d.mat', outVersion);
outFile = strrep(connFile, '.mat', outFile);

Util.saveStruct(outFile, out);
Util.protect(outFile);
