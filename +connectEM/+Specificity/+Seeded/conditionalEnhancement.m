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
curTargetId = 1;
curAxonClass = axonClasses(2);

curClassConn = classConn(curAxonClass.axonIds, :);
curHasTargetSyn = curClassConn(:, curTargetId) > 0;
curNumSyn = sum(curClassConn, 2);

%{
curHasTargetSyn = [1; 1; 0];
curNumSyn = [1; 12; 5];
%}

curPolyDegs = zeros(numel(curNumSyn), max(curNumSyn) + 1);
for curIdx = 1:size(curPolyDegs, 1)
    curNi = curNumSyn(curIdx);
    curBi = curHasTargetSyn(curIdx);
    
    if curBi
        curPolyDegs(curIdx, 1 + curNi) = 1;
    else
        curPolyDegs(curIdx, 1 + [0, curNi]) = [1, -1];
    end
end

curFun = @(x) sum((x .^ ((1:size(curPolyDegs, 2)) - 1)) .* curPolyDegs, 2);
1 - lsqnonlin(curFun, 0.5, 0, 1)
