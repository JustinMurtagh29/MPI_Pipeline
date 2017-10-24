% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

info = Util.runInfo();

% configuration
minEvidence = 54;
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
cacheDir = '/tmpscratch/amotta/l4/2017-10-23-axon-dendrite-overlap';
axonFile = fullfile(rootDir, 'chiasmataSplitting/20171009T193744-kmb-on-axons-6c/outputs/20171023T115300_results.mat');
outputFile = fullfile(rootDir, 'connectomeState', 'axons.mat');

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

axonFlights = load(fullfile(cacheDir, 'axon-flights.mat'), 'axonFlights');
axonFlights = axonFlights.axonFlights;

axons = load(axonFile);

%% find small agglomerates to pick up
axonSmallIds = find(~axons.indBigAxons);
axonAgglosSmall = axons.axons(axonSmallIds);
axonAgglosSmall = Superagglos.getSegIds(axonAgglosSmall);

maxSegId = Seg.Global.getMaxSegId(param);
axonAgglosLUT = Agglo.buildLUT(maxSegId, axonAgglosSmall);
axonAgglosLUT = cat(1, 0, axonAgglosLUT(:));

axonFlights.smallAxonId = axonAgglosLUT(1 + axonFlights.segId);
axonFlights(~axonFlights.smallAxonId, :) = [];

% pool evidence over agglomerates
[axonOverlap, ~, axonEvidence] = unique( ...
    axonFlights(:, {'axonId', 'smallAxonId'}), 'rows');
axonOverlap.evidence = accumarray(axonEvidence, 1);

% discard overlaps below evidence threshold
axonOverlap = sortrows(axonOverlap, 'evidence', 'descend');
axonOverlap(axonOverlap.evidence < minEvidence, :) = [];

% assign to axon with largest evidence
[~, uniRows] = unique(axonOverlap(:, {'axonId', 'smallAxonId'}), 'rows');
axonOverlap = axonOverlap(uniRows, :);

% build final agglomerates
axonLargeIds = find(axons.indBigAxons);
axonAgglosLarge = axons.axons(axonLargeIds);
axonAgglosLarge = Superagglos.getSegIds(axonAgglosLarge);

pickedUpAgglos = accumarray( ...
    axonOverlap.axonId, (1:size(axonOverlap, 1))', size(axonAgglosLarge), ...
    @(r) {cat(1, axonAgglosSmall{r})}, {zeros(0, 1)});

axonAgglos = cellfun( ...
    @vertcat, axonAgglosLarge, pickedUpAgglos, ...
    'UniformOutput', false);

%% save result
out = struct;
out.axonIds = axonIds;
out.axonAgglos = axonAgglos;
out.runInfo = info;

Util.saveStruct(outputFile, out);
system(sprintf('chmod a-w "%s"', outputFile));