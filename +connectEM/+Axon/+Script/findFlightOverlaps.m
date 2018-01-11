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

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

axons = load(axonFile);

%% use WKW segmentation
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% look up flight paths
if ~exist(cacheFile, 'file')
    [axonFlights, axonFlightsMeta] = ...
        Superagglos.getFlightPathSegIds( ...
            param, axons.axons(axons.indBigAxons), nhoodSize);
    save(cacheFile, '-v7.3', 'axonFlights', 'axonFlightsMeta');
else
    load(cacheFile, 'axonFlights', 'axonFlightsMeta');
end

%% what we want:
% Overlap between a pair of super-agglomerates
%
% Overlap between a pair of segment-based agglomerates should always be
% empty (by construction for segment equivalence classes).
% Overlap between flight path and segment-based agglomerates is easy
% because each segment belongs to at most one segment-based agglomerate.

%% pairwise overlap between flight-based parts of agglomerates

% remove entries with segment ID zero
axonFlights = axonFlights(axonFlights.segId > 0, :);

% remove per-node multiplicity
axonFlights = unique(axonFlights, 'rows');

% count per-agglomerate
[aggloTable, ~, uniRows] = unique( ...
    axonFlights(:, {'aggloId', 'segId'}), 'rows');
aggloTable.evidence = accumarray(uniRows, 1);

[~, ~, uniRows] = ...
    unique(aggloTable.segId);
uniSegRows = accumarray( ...
    uniRows, 1:numel(uniRows), [], ...
    @(rows) {reshape(rows, [], 1)});

% only consider segments contacted by at least two agglomerates
uniSegRows(cellfun(@numel, uniSegRows) < 2) = [];

% calculate pair-wise overlap
aggloCount = sum(axons.indBigAxons);
pairEvidence = sparse(aggloCount, aggloCount);

% looks slow, but is fast enough
for curIdx = 1:numel(uniSegRows)
    curRows = uniSegRows{curIdx};
    curAgglos = aggloTable(curRows, {'aggloId', 'evidence'});
    
    for curRowA = 1:(size(curAgglos, 1) - 1)
        for curRowB = (curRowA + 1):size(curAgglos, 1)
            curIdA = curAgglos.aggloId(curRowA);
            curIdB = curAgglos.aggloId(curRowB);
            
            % update evidence
            pairEvidence(curIdB, curIdA) = ...
                pairEvidence(curIdB, curIdA) ...
              + min(curAgglos.evidence([curRowA, curRowB]));
        end
    end
end

