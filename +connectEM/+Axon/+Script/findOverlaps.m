% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
cacheDir = '/tmpscratch/amotta/l4/2018-01-axon-16b-flight-paths';
axonFile = fullfile(rootDir, 'aggloState', 'axons_16_b.mat');

% For now, let's use the same neighborhood size as in
% `+connectEM/axonAggloStateVisualization.m`.
nhoodSize = 3;

[~, axonName] = fileparts(axonFile);
cacheFile = fullfile(cacheDir, strcat(axonName, '_flights.mat'));
clear axonName;

debugDir = fullfile(cacheDir, 'debug');
mkdir(debugDir);

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

axons = load(axonFile);
axons = axons.axons(axons.indBigAxons);

%% use WKW segmentation
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% build lookup table
numAgglos = numel(axons);
aggloSegIds = Superagglos.getSegIds(axons);
aggloLUT = Agglo.buildLUT(maxSegId, aggloSegIds);

%% look up flight paths
if ~exist(cacheFile, 'file')
   [axonFlights, axonFlightsMeta] = ...
        Superagglos.getFlightPathSegIds(param, axons, nhoodSize);
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
save(cacheFile, '-append', 'overlaps');

%% show pairs with a lot of overlap
numOverlaps = nnz(overlaps);
descPairs = nan(numOverlaps, 2);

[descPairs(:, 2), descPairs(:, 1), descOverlaps] = find(overlaps);
[descOverlaps, descRows] = sort(descOverlaps, 'descend');
descPairs = descPairs(descRows, :);

%% export examples
% To determine a good threshold we want examples from the whole range.
numBins = 1;
numPairs = 10;

binEdges = linspace(log10(20), log10(50), numBins + 1);
bins = discretize(log10(descOverlaps), binEdges);

rng(0);
samplePairs = zeros(0, 1);

for curBin = 1:numBins
    curPairIds = find(bins == curBin);
    curPairIds = curPairIds(randperm(numel(curPairIds), numPairs));
    samplePairs = cat(1, samplePairs, curPairIds(:));
end

for curIdx = 1:size(samplePairs, 1)
    curPairId = samplePairs(curIdx);
    
    curPair = descPairs(curPairId, :);
    curOverlap = descOverlaps(curPairId);
    curAgglos = axons(curPair);
    
    curSkel = skeleton();
    curSkel = curSkel.addTree( ...
        sprintf('Agglomerate %d', curPair(1)), ...
        curAgglos(1).nodes(:, 1:3), curAgglos(1).edges, [1, 0, 0, 1]);
    curSkel = curSkel.addTree( ...
        sprintf('Agglomerate %d', curPair(2)), ...
        curAgglos(2).nodes(:, 1:3), curAgglos(2).edges, [0, 0, 1, 1]);
    curSkel = curSkel.setDescription(...
        sprintf('Overlap: %d', curOverlap));
    
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    
    curSkelFile = fullfile(debugDir, sprintf( ...
        '%d_agglos-%d-and-%d.nml', curIdx, curPair));
    curSkel.write(curSkelFile);
end

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
                  + min(curAgglos.evidence([curRowA, curRowB]));
            end
        end
    end
end


