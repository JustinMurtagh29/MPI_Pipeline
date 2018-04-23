% Linearize all axons by splitting all branch points.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

chiasmaFile = fullfile(...
    rootDir, 'tripletDetection', ...
    '20180423T105503-on-axons-18b', ...
    '20180423T105912_chiasmata.mat');

%% Loading data
% Parameter
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

% Chiasmata
data = load(chiasmaFile);
chiasmata = data.chiasmata;

axonFile = data.info.param.axonFile;
chiasmaParam = data.info.param.chiasmaParam;
clear data;

% Super-agglomerates
axons = load(axonFile);
axons = axons.axons(axons.indBigAxons);
assert(isequal(size(axons), size(chiasmata)));

% NOTE(amotta): Endings and solved chiasmata were messed up in axons 18b.
% So let's just replace them with correctly shaped null values.
for curIdx = 1:numel(axons)
    axons(curIdx).endings = zeros(0, 1);
    axons(curIdx).solvedChiasma = ...
        false(size(axons(curIdx).nodes, 1), 1);
end

%% Linearize axons
aggloCount = numel(axons);
out = cell(aggloCount, 1);

tic;
for curIdx = 1:aggloCount
    curAgglos = axons(curIdx);
    curChiasmata = chiasmata{curIdx};
    
    curNrExits = curChiasmata.ccCenterIdx;
    curNrExits = curChiasmata.nrExits(curNrExits);
    
    % Fake dangling tracings
    curOverlaps = arrayfun( ...
        @(n) [transpose(1:n), zeros(n, 1)], ...
        curNrExits, 'UniformOutput', false);
    
    if ~isempty(curOverlaps)
        curAgglos = connectEM.Chiasma.Triplet.splitAgglo( ...
            param, chiasmaParam, curAgglos, curChiasmata, curOverlaps);
    end
    
    out{curIdx} = curAgglos;
    Util.progressBar(curIdx, aggloCount);
end
