% We've automatically generated axon agglomerates, queried their endings
% and applied the safe (i.e., non-E3a) flight paths. This produced the
% state axons_06_c. On these axons we then detected and solved >= 4-fold
% chiasmata. This resulted in axons_08_a.
%
% Now we try to apply the flight paths which were previously classified as
% E3a to axons_08_a. Importantly, all the overlap information contained in
% the previously calculated files refers to a different axon state!
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
minNodeEvidence = 2 * 27;
maxNumFlights = 8000;

rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
workingDir = fullfile(rootDir, 'aggloState');

usedFlightFiles = fullfile(workingDir, {'axons_10_a_plus_10kE3a_c.mat'});
flightFile = fullfile(workingDir, 'caseDistinctions_6.0.mat');
caseFile = fullfile(workingDir, 'attachedEndings_6.0.mat');
axonFile = fullfile(workingDir, 'axons_11_a.mat');

info = Util.runInfo();

%% load data
data = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = data.p;

maxSegId = Seg.Global.getMaxSegId(param);

data = load(axonFile, 'axons', 'indBigAxons');
allAxons = data.axons;
axonIds = find(data.indBigAxons);
axons = data.axons(axonIds);

data = load(flightFile, 'flightPaths');
flights = data.flightPaths.ff;

flights = structfun( ...
    @(vals) reshape(vals, [], 1), ...
    flights, 'UniformOutput', false);
flights.id = reshape(1:numel(flights.filenames), [], 1);

data = load(caseFile, 'endingCaseDistinctions', 'flightsOfEndingCases');
endingCases = data.endingCaseDistinctions;
flightsForEnding = data.flightsOfEndingCases;

%% restrict to E3a endings
endingMask = (endingCases == 5);
flightIds = cat(2, flightsForEnding{endingMask});
flightIds = unique(flightIds);

fprintf('# E3a endings: %d\n', sum(endingMask));
fprintf('# flights from E3a endings: %d\n', numel(flightIds));

flights = structfun( ...
    @(vals) vals(flightIds), ...
    flights, 'UniformOutput', false);
clear flightIds;

%% remove used flights
usedFlightIds = cellfun( ...
    @(f) load(f, 'usedFlightIds'), ...
    usedFlightFiles, 'UniformOutput', false);
usedFlightIds = Util.concatStructs('last', usedFlightIds{:});
usedFlightIds = usedFlightIds.usedFlightIds;

fprintf('\n');
fprintf('# E3a flights already used: %d\n', numel(usedFlightIds));
fprintf('⇒ These flights are being removed\n');

keepMask = cat(1, flights.id);
keepMask = not(ismember(keepMask, usedFlightIds));
clear usedFlightIds;

flights = structfun( ...
    @(vals) vals(keepMask, :), ...
    flights, 'UniformOutput', false);
clear keepMask;

%% determine overlap between flight paths and axon agglomerates
axonAgglos = Superagglos.getSegIds(axons);

flightOverlaps = ...
    connectEM.Flight.overlapWithAgglos( ...
        param, flights, axonAgglos, ...
        'minStartEvidence', 13, ...
        'minEndEvidence', 2 * 27);

axonLUT = Agglo.buildLUT(maxSegId, axonAgglos);
axonLUT = cat(1, 0, reshape(axonLUT, [], 1));

%% restrict to flights with clear overlaps
startOkay = (cellfun(@numel, flightOverlaps(:, 1)) == 1);
endOkay = (cellfun(@numel, flightOverlaps(:, 2)) == 1);
bothOkay = startOkay & endOkay;

fprintf('\n');
fprintf('# flights with clear start: %d\n', sum(startOkay));
fprintf('# flights with clear end: %d\n', sum(endOkay));
fprintf('# flights with clear start and end: %d\n', sum(bothOkay));

flights = structfun( ...
    @(vals) vals(bothOkay), ...
    flights, 'UniformOutput', false);
flights.overlaps = flightOverlaps(bothOkay, :);
flights.overlaps = cell2mat(flights.overlaps);

%% remove redundant flights
axonCount = numel(axons);

% avoid duplicate edges
uniRows = sort(flights.overlaps, 2);
[~, uniRows] = unique(uniRows, 'rows');

flights = structfun( ...
    @(vals) vals(uniRows, :), ...
    flights, 'UniformOutput', false);

%% execute only a subset of flights
rng(0);

flightCount = numel(flights.filenames);
execNumFlights = min(flightCount, maxNumFlights);

randIds = randperm(flightCount);
randIds = randIds(1:execNumFlights);
randIds = reshape(randIds, [], 1);

flights = structfun( ...
    @(vals) vals(randIds, :), ...
    flights, 'UniformOutput', false);

fprintf('\n');
fprintf('# flights selected for execution: %d\n', execNumFlights);



%% evaluation
axonStats = table;
axonStats.axonId = reshape(1:numel(out.axons), [], 1);
axonStats.nrNodes = arrayfun(@(a) size(a.nodes, 1), out.axons);
axonStats.nrSegments = arrayfun( ...
    @(a) sum(not(isnan(a.nodes(:, 4)))), out.axons);
axonStats.nrFlights = accumarray( ...
    flights.axonComp, 1, size(out.axons), @sum, 0);
axonStats = sortrows(axonStats, 'nrSegments', 'descend');

fprintf('\n');
fprintf('Largest axons:\n\n');
disp(axonStats(1:20, :));