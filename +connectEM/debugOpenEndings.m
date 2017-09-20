% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

%% configuration
% agglomeration state
% 5.0 is the one which also contains CS_MB_L4_axEndQuerySpecial2_16_09_2017
state = '5.0';

% minimum distance from end of dataset (EoD) for valid endings
minEodDistNm = 3000;

%% load parameters
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

%% build paths to data
dataDir = fullfile(param.saveFolder, 'aggloState');
[skeletonFolders, suffixFlightPaths, suffix] = connectEM.setQueryState(state);    

%% load endings and flight paths
% TODO(amotta): Find out what set of endings this is.
% Q: Were endings at the end of dataset (EoD) already filtered out?
% A: Nope, doesn't look like it.
endings = load(fullfile(dataDir, 'axonEndings.mat'));
endingOverlap = load(fullfile(dataDir, strcat('axonEndingOverlaps', suffix, '.mat')));

% TODO(amotta): Find out what set of flight paths this is.
% Q: Were NMLs with comments already filtered out?
% A: Yes. See +connectEM/getAggloQueryOverlapB.m line 23
%    Paths without valid start node are removed as well (line 24)
% Q: What about paths which reached no or multiple agglomerates?
% A: Nope, they are not removed

% TODO(amotta): I could instead use the `idxGood` field. But remember that
% this does not yet take care of multiple ending agglomerates
flightPaths = load(fullfile(dataDir, strcat('axonPostQueryAnalysisState', suffix, '.mat')), 'ff', 'startAgglo', 'endAgglo');

%% count all endings
endingCount = sum(cellfun(@max, endings.borderClusters));
fprintf('Number of endings: %d\n', endingCount);

%% count endings in core
% first, determine ending positions (as center of mass of borders)
borderIdsOff = cumsum(cellfun(@max, endings.borderClusters));
borderIdsOff = num2cell(cat(1, 0, borderIdsOff(1:(end - 1))));

borderIds = cellfun( ...
    @(ids, off) ids + off, ...
    endings.borderClusters, borderIdsOff, ...
    'UniformOutput', false);
borderIds = cell2mat(borderIds);
borderPositions = cell2mat(endings.borderPositions);

endingPositions = accumarray( ...
    borderIds, 1:numel(borderIds), ...
    [], @(rows) {mean(borderPositions(rows, :), 1)});
endingPositions = cell2mat(endingPositions);

% next, check which of these are within core
minEodOff = ceil(minEodDistNm ./ param.raw.voxelSize);

boxEod = param.bbox;
boxEod = cat(2, ...
    boxEod(:, 1) + minEodOff(:), ...
    boxEod(:, 2) - minEodOff(:));

% then find endings in core
endingsInCore = ...
    all(bsxfun(@ge, endingPositions, boxEod(:, 1)'), 2) ...
  & all(bsxfun(@le, endingPositions, boxEod(:, 2)'), 2);
endingsInCore = find(endingsInCore);
endingCount = numel(endingsInCore);

fprintf('Number of endings within dataset core: %d\n', endingCount);

%% count seeded endings
seededEndings = endingOverlap.startEndingOverlaps;
seededEndings = cat(1, seededEndings{:});

% make unique and remove invalid endings,
% i.e., no-ending and end-of-dataset endings
seededEndings = intersect(seededEndings, endingsInCore);
seededEndingCount = numel(seededEndings);

fprintf('\n');
fprintf('Number of seeded endings: %d\n', seededEndingCount);
fprintf('Number of non-seeded endings: %d\n', endingCount - seededEndingCount);

%% cout reached endings
reachedEndings = endingOverlap.endEndingOverlaps;
reachedEndings = cat(1, reachedEndings{:});

% make unique and remove zero (i.e., no ending)
reachedEndings = intersect(reachedEndings, endingsInCore);
reachedEndingCount = numel(reachedEndings);

fprintf('\n');
fprintf('Number of reached endings: %d\n', reachedEndingCount);
fprintf('Number of non-reached endings: %d\n', endingCount - reachedEndingCount);

%% count encountered endings (be it as seed or end)
encounteredEndings = cat(1, seededEndings, reachedEndings);
encounteredEndings = unique(encounteredEndings);
encounteredEndingCount = numel(encounteredEndings);

fprintf('\n');
fprintf('Number of encountered endings: %d\n', encounteredEndingCount);
fprintf('Number of non-encountered endings: %d\n', endingCount - encounteredEndingCount);