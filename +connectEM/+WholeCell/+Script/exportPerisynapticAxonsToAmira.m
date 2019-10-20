% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outDir = '/tmpscratch/amotta/l4/2019-10-20-whole-cell-input-axon-snippet-isosurfaces';

periSynRadNm = 10E3;

isoParams = { ...
    'reduce', 0.05, ...
    'smoothSizeHalf', 4, ...
    'smoothWidth', 8};

segParam = struct;
segParam.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
segParam.backend = 'wkwrap';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.seg = segParam;

segPos = Seg.Global.getSegToPointMap(param);
[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Find synapses onto whole cells
clear cur*;

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.cellId = conn.denMeta.cellId(synT.postAggloId);

synT = synT(synT.cellId > 0, :);
synT = sortrows(synT, {'cellId', 'id'});

%% Build axon snippet for each synapse
clear cur*;

curBoxPad = ceil(periSynRadNm ./ param.raw.voxelSize(:));

tic;
synT.axonSnippet(:) = {[]};
for curIdx = 1:height(synT)
    curSynSegIds = syn.synapses.presynId{synT.id(curIdx)};
    curSynBox = segPos(curSynSegIds, :);
    curSynBox = [ ...
        transpose(min(curSynBox, [], 1)), ...
        transpose(max(curSynBox, [], 1))] ...
        + [-1, +1] .* curBoxPad;

    curAxonId = synT.preAggloId(curIdx);
    curAxonSegIds = conn.axons{curAxonId};

    curMask = segPos(curAxonSegIds, :);
    curMask = ...
        all(curSynBox(:, 1)' <= curMask, 2) ...
      & all(curMask <= curSynBox(:, 2)', 2);
    curAxonSegIds = curAxonSegIds(curMask);

    [~, curSynIds] = ismember(curSynSegIds, curAxonSegIds);
    curSynIds = setdiff(curSynIds, 0);

    curAxon = segPos(curAxonSegIds, :) .* param.raw.voxelSize;
    curAxon = squareform(pdist(curAxon));
    curAxon(1, 1) = 0;

    curAxon = graph(curAxon);
    curAxon = minspantree(curAxon);
    curMask = distances(curAxon, curSynIds);
    curMask = min(curMask, [], 1) <= periSynRadNm;

    curAxonSegIds = curAxonSegIds(curMask);
    synT.axonSnippet{curIdx} = curAxonSegIds(:);

    Util.progressBar(curIdx, height(synT));
end

%% Save results
clear cur*;

curOut = struct;
curOut.info = info;
curOut.synT = synT;

curOutFile = fullfile(outDir, 'input-synapses.mat');
Util.saveStruct(curOutFile, curOut);
Util.protect(curOutFile);

%% Build isosurfaces
clear cur*;

param.seg = segParam;
Visualization.exportAggloToAmira( ...
    param, synT.axonSnippet, outDir, isoParams{:});
