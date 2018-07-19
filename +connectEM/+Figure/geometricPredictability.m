% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
availFile = '/tmpscratch/amotta/l4/2018-07-18-surface-availability-for-connectome-v7-partially-split/axon-availability_v2.mat';

targetClasses = { ...
    'Somata', 'SO';
    'ProximalDendrite', 'PD'; ...
    'SmoothDendrite', 'SD'; ...
    'ApicalDendrite', 'AD'; ...
    'AxonInitialSegment', 'AIS'; ...
    'OtherDendrite', 'Other'};

targetTags = reshape(targetClasses(:, 2), 1, []);
targetClasses = reshape(targetClasses(:, 1), 1, []);

minSynPre = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

avail = load(availFile);

% Remove "all" class (what a stupid idea that was...)
avail.axonAvail(avail.targetClasses == 'All', :, :) = [];
avail.targetClasses(avail.targetClasses == 'All') = [];

%% Prepare data
% Class connectome
[classConn, classIds] = ...
    connectEM.Connectome.buildClassConnectome(conn);

% Fix order
assert(numel(classIds) == numel(targetClasses));
[~, classIds] = ismember(targetClasses, classIds);
classConn = classConn(:, classIds);

% Availabilities
availabilities = avail.axonAvail;
availabilities = availabilities ./ sum(availabilities, 1);

assert(numel(avail.targetClasses) == numel(targetClasses));
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = availabilities(classIds, :, :);

%% Prototype
clear cur*;

curAxonIds = axonClasses(2).axonIds;
curTargetId = 2;
curDistId = 10;

curClassConn = classConn(curAxonIds, :);
curSynCountAll = sum(curClassConn, 2);
curSynTargetAll = curClassConn(:, curTargetId);
curSynFracs = curSynTargetAll ./ curSynCountAll;

curAvails = availabilities(:, curDistId, curAxonIds);
curAvails = reshape(curAvails, size(curAvails, 1), []);
curAvails = transpose(curAvails);

curS = size(targetClasses);
curOpts = optimoptions('lsqlin', 'Display', 'off');

% Find optimal parameters
[curParam, ~, ~, curExitFlag] = lsqlin( ...
    curAvails, curSynFracs, ones(curS), 1, [], [], ...
    zeros(curS), ones(curS), ones(curS) / numel(curS), curOpts);
assert(curExitFlag == 1);
