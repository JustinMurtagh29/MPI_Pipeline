% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

minSynPre = 10;

targetClasses = { ...
    'Somata', 'ProximalDendrite', 'ApicalDendrite', ...
    'SmoothDendrite', 'AxonInitialSegment', 'OtherDendrite'};

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);

[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

classConn = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', targetClasses);

%% Testing
clear cur*;
curAxonClass = axonClasses(2);
curClassConn = classConn(curAxonClass.axonIds, :);

curSyns = sum(curClassConn, 2);

for curTargetId = 1:numel(targetClasses)
    curHit = curClassConn(:, curTargetId) > 0;

    curFuncs = @(x) ...
        (1 - curHit) + (2 .* curHit - 1) .* (x .^ curSyns);
    curDerivs = @(x) ...
        (2 .* curHit - 1) .* curSyns .* (x .^ (curSyns - 1));
    curOptions = optimoptions( ...
        'lsqnonlin', 'SpecifyObjectiveGradient', true, 'display', 'off');
    curPInv = lsqnonlin(@(x) ...
        deal(curFuncs(x), curDerivs(x)), 0.5, 0, 1, curOptions);
    curP = 1 - curPInv;
    
    curBulkP = sum(curClassConn(:, curTargetId)) / sum(curClassConn(:));
    
    [curP, curBulkP] %#ok
end
