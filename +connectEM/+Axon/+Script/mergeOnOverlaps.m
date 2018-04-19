% This script looks for pair-wise overlaps between super-agglomerates. By
% definition, the segment-based parts of the agglomerates are pair-wise
% disjoint.
%
% However, the flight path part of a super-agglomerate might overlap with
% a) the segment-based part, or b) the flight path-based part of another
% super-agglomerate.
%
% To identify these overlaps, this script first maps the flight nodes to
% segment IDs. It then separately calculates the pair-wise a) flight-to-
% flight, and b) flight-to-segment overlaps.
%
% These two contributions are then added about to obtain a single overlap
% evidence value per pair of super-agglomerates.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
cacheDir = '/tmpscratch/amotta/l4/2018-01-19-axon-17a-flight-paths';
axonFile = fullfile(rootDir, 'aggloState', 'axons_17_a.mat');
outFile = fullfile(rootDir, 'aggloState', 'axons_18_b.mat');

% For now, let's use the same neighborhood size as in
% `+connectEM/axonAggloStateVisualization.m`.
nhoodSize = 3;

[~, axonName] = fileparts(axonFile);
cacheFile = fullfile(cacheDir, strcat(axonName, '_flights.mat'));
clear axonName;

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

temp = load(axonFile);
allAxons = temp.axons(:);
largeAxonIds = find(temp.indBigAxons(:));
smallAxonIds = find(~temp.indBigAxons(:));
clear temp;

% sort edges
% I've messed up in axons 16b.
for curIdx = 1:numel(allAxons)
    allAxons(curIdx).edges = ...
        sort(allAxons(curIdx).edges, 2);
end

%% use WKW segmentation
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% build lookup table
numAgglos = numel(largeAxonIds);
aggloSegIds = Superagglos.getSegIds(allAxons(largeAxonIds));
aggloLUT = Agglo.buildLUT(maxSegId, aggloSegIds);

%% look up flight paths
if ~exist(cacheFile, 'file')
   [axonFlights, axonFlightsMeta] = ...
        Superagglos.getFlightPathSegIds( ...
            param, allAxons(largeAxonIds), nhoodSize);
    save(cacheFile, '-v7.3', 'axonFlights', 'axonFlightsMeta');
else
    load(cacheFile, 'axonFlights', 'axonFlightsMeta');
end

%% simplify / process flight table
% remove entries with segment ID zero
axonFlights(~axonFlights.segId, :) = [];

% remove per-node multiplicity
axonFlights = unique(axonFlights, 'rows');

% count per-agglomerate
[aggloTable, ~, uniRows] = unique( ...
    axonFlights(:, {'aggloId', 'segId'}), 'rows');
aggloTable.evidence = accumarray(uniRows, 1);

%% calculate overlaps
faOverlaps = flightAggloOverlaps(numAgglos, aggloLUT, aggloTable);
ffOverlaps = flightFlightOverlaps(numAgglos, aggloTable);

overlaps = faOverlaps + ffOverlaps;

%% show pairs with a lot of overlap
numOverlaps = nnz(overlaps);
descPairs = nan(numOverlaps, 2);

[descPairs(:, 2), descPairs(:, 1), descOverlaps] = find(overlaps);
[descOverlaps, descRows] = sort(descOverlaps, 'descend');
descPairs = descPairs(descRows, :);

%% find axons to be merged
execOverlaps = (overlaps >= 25);
[compCount, compLUT] = graphconncomp( ...
    execOverlaps, 'Directed', false);
compLUT = reshape(compLUT, [], 1);

tic;
disp('Merging axons...');

out = struct;
out.axons = allAxons([]);

for curIdx = 1:compCount
    curAxonIds = largeAxonIds(compLUT == curIdx);
    curAxons = allAxons(curAxonIds);
    
    if isscalar(curAxons)
        out.axons(curIdx) = curAxons;
        continue;
    end
    
    curAxonLens = Superagglos.mstLength( ...
        curAxons, param.raw.voxelSize);
   [~, curSeedAxonId] = min(curAxonLens);
    
    % find minimal set of edges to execute
    curOverlaps = execOverlaps(curAxonIds, curAxonIds);
    curOverlaps = graphminspantree(curOverlaps, curSeedAxonId);
    
    curEdges = nan(nnz(curOverlaps), 2);
   [curEdges(:, 2), curEdges(:, 1)] = find(curOverlaps);
   
    % sort edges by ascending size of the larger axon
    curSortIds = reshape(curAxonLens(curEdges), [], 2);
   [~, curSortIds] = sort(max(curSortIds, [], 2), 'ascend');
    curEdges = curEdges(curSortIds, :);
    
    curAxon = curAxons(curSeedAxonId);
    curMerged = false(size(curAxons));
    curMerged(curSeedAxonId) = true;
    
    while ~all(curMerged)
        curAggloIdx = find(xor( ...
            curMerged(curEdges(:, 1)), ...
            curMerged(curEdges(:, 2))), 1);
        
        curAggloIdx = curEdges(curAggloIdx, :);
        curAggloIdx(curMerged(curAggloIdx)) = [];
        assert(isscalar(curAggloIdx));
        
        curAxon = Superagglos.mergeOnOverlaps( ...
            curAxons(curAggloIdx), curAxon, ...
            'scale', param.raw.voxelSize, ...
            'overlapDistNm', 200, ...
            'minLenNm', 2000);
        
        curMerged(curAggloIdx) = true;
    end
    
    out.axons(curIdx) = curAxon;
    Util.progressBar(curIdx, compCount);
end
toc;

out.axons = vertcat( ...
    out.axons(:), allAxons(smallAxonIds));
SuperAgglo.check(out.axons);

out.indBigAxons = false(numel(out.axons), 1);
out.indBigAxons(1:compCount) = true;

out.parentIds = nan(numel(allAxons), 1);
out.parentIds(largeAxonIds) = compLUT;
out.parentIds(smallAxonIds) = ...
    (compCount - 1) + (1:numel(smallAxonIds));

out.info = info;

Util.saveStruct(outFile, out);
Util.protect(outFile);
    
%% overlap between flights and segment-based agglomerates
function overlap = flightAggloOverlaps(numAgglos, aggloLUT, aggloTable)
    aggloTable.aggloId(:, 2) = aggloLUT(aggloTable.segId);
    aggloTable.aggloId = sort(aggloTable.aggloId, 2);
    
    aggloTable(~aggloTable.aggloId(:, 1), :) = [];
    aggloTable(~diff(aggloTable.aggloId, 1, 2), :) = [];
    
   [~, uniRows, uniIds] = unique(aggloTable.aggloId, 'rows');
   
    pairTable = aggloTable(uniRows, :);
    pairTable.evidence = accumarray(uniIds, aggloTable.evidence);
    
    overlap = sparse( ...
        pairTable.aggloId(:, 2), ...
        pairTable.aggloId(:, 1), ...
        pairTable.evidence, ...
        numAgglos, numAgglos);
end

%% pairwise overlap between flight-based parts of agglomerates
function overlap = flightFlightOverlaps(numAgglos, aggloTable)
   [~, ~, uniRows] = ...
        unique(aggloTable.segId);
    uniSegRows = accumarray( ...
        uniRows, 1:numel(uniRows), [], ...
        @(rows) {reshape(rows, [], 1)});

    % only consider segments contacted by at least two agglomerates
    uniSegRows(cellfun(@numel, uniSegRows) < 2) = [];

    % calculate pair-wise overlap
    overlap = sparse(numAgglos, numAgglos);

    % looks slow, but is fast enough
    for curIdx = 1:numel(uniSegRows)
        curRows = uniSegRows{curIdx};
        curAgglos = aggloTable(curRows, {'aggloId', 'evidence'});

        for curRowA = 1:(size(curAgglos, 1) - 1)
            for curRowB = (curRowA + 1):size(curAgglos, 1)
                curIdA = curAgglos.aggloId(curRowA);
                curIdB = curAgglos.aggloId(curRowB);
                
                % `curIdB > curIdA` because `curRowB > curRowA`, `curRows`
                % is monotonically increasing and `aggloTable` is sorted in
                % ascending order. Hence we're only working with the
                % lower-diagonal part of `pairEvidence`.

                % update evidence
                overlap(curIdB, curIdA) = ...
                    overlap(curIdB, curIdA) ...
                  + min(curAgglos.evidence([curRowA, curRowB])); %#ok
            end
        end
    end
end
