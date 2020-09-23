% Linearize all axons by splitting all branch points.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

chiasmaFile = fullfile( ...
    rootDir, 'tripletDetection', ...
    '20200608T182413-on-axons-04', ...
    '20200608T182615_chiasmata.mat');

info = Util.runInfo();
Util.showRunInfo(info);

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

% Output file
[outDir, outFile] = fileparts(axonFile);
outFile = fullfile(outDir, sprintf('%s_linearized.mat', outFile));

% Super-agglomerates
axons = load(axonFile);

% NOTE(amotta): This early in the axon reconstruction, we haven't set the
% `endings` and `solvedChiasma` fields yet. So, let's just set them to
% correctly shaped null values.
for curIdx = 1:numel(axons.axons)
    axons.axons(curIdx).endings = zeros(0, 1);
    axons.axons(curIdx).solvedChiasma = ...
        false(size(axons.axons(curIdx).nodes, 1), 1);
end

bigAxonIds = find(axons.indBigAxons);
bigAxons = axons.axons(bigAxonIds);

% Sanity check
assert(isequal(size(bigAxons), size(chiasmata)));

%% Linearize axons
aggloCount = numel(bigAxons);
split = cell(aggloCount, 1);

tic;
for curIdx = 1:aggloCount
    curAgglos = bigAxons(curIdx);
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
    
    split{curIdx} = curAgglos;
    Util.progressBar(curIdx, aggloCount);
end

%% Build output structure
out = struct;

% Big axons
out.axons = cat(1, split{:});
out.parentIds = repelem(bigAxonIds, cellfun(@numel, split));

% Small axons
smallAxonIds = find(~axons.indBigAxons);
out.axons = [out.axons; axons.axons(smallAxonIds)];
out.parentIds = [out.parentIds; smallAxonIds];

% Completing
out.indBigAxons = axons.indBigAxons(out.parentIds);
out.info = info;

%% Save result
Util.saveStruct(outFile, out);
Util.protect(outFile);
