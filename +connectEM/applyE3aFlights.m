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
maxNumFlights = inf;

rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
workingDir = fullfile(rootDir, 'aggloState');

axonFile = fullfile(workingDir, 'axons_10_a.mat');
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
    curEndId = axonLUT(curSegIds);
    
    % find start agglo (> 1/2 node evidence)
   [curStartAxon, ~, curStartEvidence] = unique(curEndId(1, :));
    curStartAxon = curStartAxon(accumarray(curStartEvidence(:), 1) > 13);
    curStartAxon = setdiff(curStartAxon, 0);
    
    % find end agglos
   [curEndAxons, ~, curEndEvidence] = unique(curEndId(:));
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

adjEdges = flights.overlaps;
adjEdges = sort(adjEdges, 2);

% avoid duplicate edges
[~, uniRows] = unique(adjEdges, 'rows');

adjEdges = adjEdges(uniRows, :);
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

%% grouping agglomerates
adjMat = sparse( ...
    adjEdges(:, 2), adjEdges(:, 1), ...
    true, numel(axons), numel(axons));
[axonCompCount, axonComps] = ...
    graphconncomp(adjMat, 'Directed', false);

axonComps = reshape(axonComps, [], 1);
flights.axonComp = axonComps(flights.overlaps(:, 1));

%% patch in flight paths
out = struct;
out.info = info;
out.axons = axons([]);

for curComp = 1:axonCompCount
    curAxons = (axonComps == curComp);
    curAxons = axons(curAxons);
    
    if numel(curAxons) == 1
        out.axons(curComp) = curAxon;
        continue;
    end
    
    curAxon = struct;
    curAxon.nodes = cat(1, curAxons.nodes);
    
    curSegLUT = zeros(maxSegId, 1);
    curNodeIds = find(not(isnan(curAxon.nodes(:, 4))));
    curSegLUT(curAxon.nodes(curNodeIds, 4)) = curNodeIds;
    
    % determine node offset
    curNodeOff = arrayfun(@(a) size(a.nodes, 1), curAxons);
    curNodeOff = cumsum(cat(1, 0, curNodeOff(1:(end - 1))));
    
    curAxon.edges = cell2mat(arrayfun(@(ax, off) ...
        ax.edges + off, curAxons, curNodeOff, 'UniformOutput', false));
    curAxon.endings = cell2mat(arrayfun(@(ax, off) ...
        ax.endings + off, curAxons, curNodeOff, 'UniformOutput', false));
    
    % patch in flight paths
    curFlightIds = find(flights.axonComp == curComp);
    curFlightIds = reshape(curFlightIds, 1, []);
    
    curFlightNodes = cell(numel(curFlightIds), 1);
    curFlightEdges = cell(numel(curFlightIds), 1);
    curNodeOff = size(curAxon.nodes, 1);
    
    for curIdx = 1:numel(curFlightIds)
        curId = curFlightIds(curIdx);
        
        curSegIds = cat(2, ...
            flights.segIds{curId}, ...
            flights.neighbours{curId});
        curSegIds = transpose(curSegIds);
        
        curEndId = axonLUT(1 + curSegIds);
        curOverlaps = flights.overlaps(curId, :);
        
        % find flight path stretch to extract
        curStartMask = any(curEndId == curOverlaps(1), 1);
        curEndMask = any(curEndId == curOverlaps(2), 1);
        
        curEndIdx = 1 + find(curEndMask(2:end), 1, 'first');
        curStartIdx = max(find(curStartMask(1:(curEndIdx - 1)))); %#ok
        
        curStartNodeId = find( ...
            curEndId(:, curStartIdx) == curOverlaps(1), 1);
        curStartNodeId = curSegIds(curStartNodeId, curStartIdx);
        curStartNodeId = curSegLUT(curStartNodeId);
        
        curEndNodeId = find( ...
            curEndId(:, curEndIdx) == curOverlaps(2), 1);
        curEndNodeId = curSegIds(curEndNodeId, curEndIdx);
        curEndNodeId = curSegLUT(curEndNodeId);
        
        % path in flight nodes
        curNodesToAdd = flights.nodes{curId};
        curNodesToAdd(:, 4) = nan;
        
        % truncate flights
        curNodesToAdd = curNodesToAdd( ...
            (curStartIdx + 1):(curEndIdx - 1), :);
        
        % collect new edges
        curNodeCount = size(curNodesToAdd, 1);
        curEdgesToAdd = zeros((curNodeCount - 1) + 2, 2);
        curEdgesToAdd((1 + 1):end, 1) = 1:curNodeCount;
        curEdgesToAdd(1:(end - 1), 2) = 1:curNodeCount;
        
        curEdgesToAdd = curEdgesToAdd + curNodeOff;
        curNodeOff = curNodeOff + curNodeCount;
        
        % connect to super-agglomerates
        curEdgesToAdd(1) = curStartNodeId;
        curEdgesToAdd(end) = curEndNodeId;
        
        curFlightNodes{curIdx} = curNodesToAdd;
        curFlightEdges{curIdx} = curEdgesToAdd;
    end
    
    curAxon.nodes = cat(1, curAxon.nodes, cell2mat(curFlightNodes));
    curAxon.edges = cat(1, curAxon.edges, cell2mat(curFlightEdges));
    curAxon.edges = sort(curAxon.edges, 2);
    
    curFlightNodeCount = sum(cellfun( ...
        @(n) size(n, 1), curFlightNodes));
    curAxon.solvedChiasma = cat( ...
        1, curAxons.solvedChiasma, ...
        false(curFlightNodeCount, 1));
    
    % sanity check
    curAdj = sparse( ...
        curAxon.edges(:, 2), curAxon.edges(:, 1), ...
        true, size(curAxon.nodes, 1),size(curAxon.nodes, 1));
    assert(graphconncomp(curAdj, 'Directed', false) == 1);
    
    out.axons(curComp) = curAxon;
end

otherAxonIds = setdiff(1:numel(allAxons), axonIds);
out.axons = cat(1, out.axons(:), allAxons(otherAxonIds));

out.indBigAxons = false(size(out.axons));
out.indBigAxons(1:axonCompCount) = true;

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