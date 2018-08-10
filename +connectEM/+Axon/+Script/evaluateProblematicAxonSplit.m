% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

minSynPre = 10;
info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

% NOTE(amotta): The version (v2) that was used to generate the connectome
% did not track the ancestry of axons. Here, we're using an updated version
% (v4) that contains the exact same axon super-agglomerates, but some
% additional meta data that allows for the following analyses.
axonFile = conn.info.param.axonFile;
axonFile = strrep(axonFile, '_v2.mat', '_v4.mat');
axons = load(axonFile);

%% Evaluate
numOldAxons = numel(unique(axons.parentIds(axons.parentIdsSplit ~= 0))) %#ok
numNewAxons = sum(axons.parentIdsSplit ~= 0) %#ok

splitMask = axons.parentIdsSplit(conn.axonMeta.parentId) ~= 0;
numNewAxonsInConn = sum(splitMask) %#ok

splitSynMask = splitMask & conn.axonMeta.synCount >= 10;
numNewAxonsInConnWithTenSynapses = sum(splitSynMask) %#ok

axonClassNames = categories(conn.axonMeta.axonClass);
axonClassNames = reshape(axonClassNames, 1, []) %#ok
numNewAxonsPerClass = transpose(accumarray( ...
    double(conn.axonMeta.axonClass(splitSynMask)), ...
    1, [numel(axonClassNames), 1])) %#ok
