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

rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
workingDir = fullfile(rootDir, 'aggloState');

axonFile = fullfile(workingDir, 'axons_08_a.mat');
flightFile = fullfile(workingDir, 'caseDistinctions_6.0.mat');
caseFile = fullfile(workingDir, 'attachedEndings_6.0.mat');

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

%% determine overlap between flight paths and axon agglomerates
axonAgglos = Superagglos.getSegIds(axons);
axonLUT = Agglo.buildLUT(maxSegId, axonAgglos);
axonLUT = vertcat(0, axonLUT(:));

flightCount = numel(flights.filenames);
flightOverlaps = cell(flightCount, 2);

for curIdx = 1:flightCount
    curSegIds = 1 + horzcat( ...
        flights.segIds{curIdx}, ...
        flights.neighbours{curIdx});
    curAxonIds = axonLUT(curSegIds);
    
    % find start agglo (> 1/2 node evidence)
   [curStartAxon, ~, curStartEvidence] = unique(curAxonIds(1, :));
    curStartAxon = curStartAxon(accumarray(curStartEvidence(:), 1) > 13);
    curStartAxon = setdiff(curStartAxon, 0);
    
    % find end agglos
   [curEndAxons, ~, curEndEvidence] = unique(curAxonIds(:));
    curEndAxons = curEndAxons(accumarray( ...
        curEndEvidence(:), 1) >= minNodeEvidence);
    curEndAxons = setdiff(curEndAxons, cat(1, 0, curStartAxon(:)));
    
    flightOverlaps{curIdx, 1} = curStartAxon;
    flightOverlaps{curIdx, 2} = curEndAxons;
end

%% restrict to flights with clear overlaps
startOkay = (cellfun(@numel, flightOverlaps(:, 1)) == 1);
endOkay = (cellfun(@numel, flightOverlaps(:, 2)) == 1);
bothOkay = startOkay & endOkay;

fprintf('# flights with clear start: %d\n', sum(startOkay));
fprintf('# flights with clear end: %d\n', sum(endOkay));
fprintf('# flights with clear start and end: %d\n', sum(bothOkay));

flights = structfun( ...
    @(vals) vals(bothOkay), ...
    flights, 'UniformOutput', false);

%%