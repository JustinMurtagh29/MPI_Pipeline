% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

%% configuration
param = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat', 'p');
param = param.p;
state = '6.2';

%% preamble copied from connectEM.createNewSuperagglos
dataDir = fullfile(param.saveFolder, 'aggloState');

[~, ~, suffix, axonVersion, axonVersionNew, casesToMerge] = connectEM.setQueryState(state);

% Load current state of agglomerates
agglos = load(fullfile(dataDir, strcat('axons_',num2str(axonVersion,'%.2i'),'.mat')));
superAgglos = agglos.axons(agglos.indBigAxons);

% NOTE(amotta): `agglos` contains all super-agglomerates. `superAgglos`
% is restricted to the large-enough ones.

% Load linkages and cases for execusion
load(fullfile(dataDir, strcat('caseDistinctions',suffix,'.mat')),...
    'linkagesAgglos','flightPaths','linkagesFlat');
load(fullfile(dataDir, strcat('attachedEndings',suffix,'.mat')),'flightsOfEndingCases','endingCaseDistinctions');

% Choose cases for merging and generation of queries
if isempty(casesToMerge)
    casesToMerge = [1:6, 8:14];
end
executedFlightPaths = flightsOfEndingCases(ismember(endingCaseDistinctions, casesToMerge));
executedFlightPaths = unique(cat(2,executedFlightPaths{:})');

% sanity checks
assert(isequal(size(linkagesAgglos), size(linkagesFlat)));
assert(size(linkagesAgglos, 1) == numel(flightPaths.startAgglo));
assert(size(linkagesAgglos, 1) == numel(flightPaths.endAgglo));
assert(size(linkagesAgglos, 1) == numel(flightPaths.ff.segIds));
assert(~any(diff(structfun(@numel, flightPaths.ff))));
assert(max(executedFlightPaths) < size(linkagesAgglos, 1));

%% redundancy removal
mask = false(size(linkagesAgglos, 1), 1);
mask(executedFlightPaths) = true;

% mask for flight that attach at both ends
validAttach = @(ids) numel(ids) == 1 && ids(1) > 0;
aggloPairMask = ...
    cellfun(validAttach, flightPaths.startAgglo) ...
  & cellfun(validAttach, flightPaths.endAgglo);

dangIds = find(mask & ~aggloPairMask);
aggloPairIds = find(mask & aggloPairMask);

aggloPairs = horzcat( ...
    cell2mat(flightPaths.startAgglo(aggloPairIds)), ...
    cell2mat(flightPaths.endAgglo(aggloPairIds)));
aggloPairs = sort(aggloPairs, 2);

% NOTE(amotta): Build a table with the pairs-of-agglomerates and
% pairs-of-endings. By sorting this table DESCENDING ENDING IDS for each
% pair of agglomerates and then taking the unique agglomerate pairs, we
% favor ending attachments.
endPairs = linkagesFlat(aggloPairIds, :);
endPairs = sort(endPairs, 2, 'descend');

aggloEndPairs = horzcat(aggloPairs, endPairs);
[aggloEndPairs, sortIds] = sortrows(aggloEndPairs, [1, 2, -3, -4]);
aggloPairIds = aggloPairIds(sortIds);

[uniAggloPairs, pairFlightIds, pairFlightGroups] = ...
    unique(aggloEndPairs(:, 1:2), 'rows', 'first');
pairFlightIds = aggloPairIds(pairFlightIds);

% endings for dangling ones
assert(all(linkagesFlat(dangIds, 1) > 0));
assert(all(isnan(linkagesFlat(dangIds, 2))));
dangEndIds = linkagesFlat(dangIds, 1);

[~, dangFlightIds, dangFlightGroups] = unique(dangEndIds);
dangFlightIds = dangIds(dangFlightIds);

% put it all together
uniFlightIds = sort(cat(1, pairFlightIds, dangFlightIds));
eqNew = Graph.findConnectedComponents(uniAggloPairs, true, true);
eqNew = [eqNew; num2cell(setdiff( ...
    1 : length(superAgglos), cell2mat(eqNew)))'];
