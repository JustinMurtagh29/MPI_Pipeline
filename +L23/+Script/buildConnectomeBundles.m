% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191203T021242-results_20191203T021242-results-auto-spines-v2_SynapseAgglomerates--20191203T021242-results--20191203T021242-results-auto-spines-v2--v1.mat');

outSuffix = '_with-synapses_v1';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = fullfile(rootDir, 'allParameter.mat');
param = Util.load(param, 'p');

maxSegId = Seg.Global.getMaxSegId(param);

conn = load(connFile);
syn = load(conn.info.param.synFile);

axonFile = conn.info.param.axonFile;
axons = Util.load(axonFile, 'axons');

dendriteFile = conn.info.param.dendriteFile;
dendrites = Util.load(dendriteFile, 'dendrites');

%% Calculate synapse positions
clear cur*;

curCalcPos = @(where) ...
    Synapse.calculatePositions( ...
        param, syn.synapses, where); 

curPosPre = curCalcPos('preRecenter');
curPosSyn = curCalcPos('border');
curPosPost = curCalcPos('postRecenter');

synPos = cat(2, ...
    reshape(transpose(curPosPre), 3, 1, []), ...
    reshape(transpose(curPosSyn), 3, 1, []), ...
    reshape(transpose(curPosPost), 3, 1, []));
synPos = uint32(round(synPos));

%% Build axon bundle
clear cur*;

curOut = SuperAgglo.toBundle(maxSegId, axons);
[curSynOff, curSynPos] = buildSynapses(conn, synPos, 'pre');
curOut.agglomerate_to_synapses_offsets = curSynOff;
curOut.agglomerate_to_synapses_positions = curSynPos;

[curOutDir, curOutFile] = fileparts(axonFile);
curOutFile = sprintf('%s_axon-bundle%s.hdf5', curOutFile, outSuffix);
curOutFile = fullfile(curOutDir, curOutFile);

structToHdf5(curOutFile, '/', curOut, false);
infoToHdf5(curOutFile, info);
Util.protect(curOutFile);

%% Build dendrite bundle
clear cur*;

curOut = SuperAgglo.toBundle(maxSegId, dendrites);
[curSynOff, curSynPos] = buildSynapses(conn, synPos, 'post');
curOut.agglomerate_to_synapses_offsets = curSynOff;
curOut.agglomerate_to_synapses_positions = curSynPos;

[curOutDir, curOutFile] = fileparts(dendriteFile);
curOutFile = sprintf('%s_dendrite-bundle%s.hdf5', curOutFile, outSuffix);
curOutFile = fullfile(curOutDir, curOutFile);

structToHdf5(curOutFile, '/', curOut, false);
infoToHdf5(curOutFile, info);
Util.protect(curOutFile);

%% Utilities
function [off, pos] = buildSynapses(conn, synPos, preOrPost)
    preOrPost = lower(preOrPost);
    
    switch preOrPost
        case 'pre'
            column = 1;
            count = numel(conn.axons);
        case 'post'
            column = 2;
            count = numel(conn.dendrites);
        otherwise
            error('Not reachable');
    end
            
    conn = conn.connectome;
    targetIds = cellfun(@numel, conn.synIdx);
    targetIds = repelem(conn.edges(:, column), targetIds);
    synIds = cell2mat(conn.synIdx);

    synIds = accumarray( ...
        targetIds, synIds, [count, 1], ...
        @(ids) {ids(:)}, {zeros(0, 1, 'like', synIds)});

    off = cellfun(@numel, synIds);
    off = uint64(cat(1, [0; 0], cumsum(off(:))));

    synIds = cell2mat(synIds);
    pos = uint32(synPos(:, :, synIds));
end
