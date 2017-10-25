% This function takes a bunch of axon super-agglomerates and checks if the
% flight paths contained therein have significant overlaps with some of the
% dendrites agglomerates.
%
% The hope is that this helps us to identify axon fragments which ended up
% in the wrong bucket. It could be that these fragments are a major source
% of "dendrite" endings. Identifying these fragments could potentially help
% us to get rid of unnecessary queries.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear();

% configuration
minEvidence = 10 * 27;
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'chiasmataSplitting/20171009T193744-kmb-on-axons-6c/outputs/20171023T115300_results.mat');
dendFile = fullfile(rootDir, 'aggloState/dendrites_03_v2_splitmerged.mat');

%% load input data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segSizes = Seg.Global.getSegToSizeMap(param);

axons = load(axonFile);
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

dends = load(dendFile);
dendIds = find(dends.indBigDends);
dends = dends.dendrites(dendIds);

%% prepare for dendrite overlap
dendSegIds = Superagglos.getSegIds(dends);
dendVols = cellfun(@(ids) max([0; sum(segSizes(ids))]), dendSegIds);

dendLUT = Agglo.buildLUT(maxSegId, dendSegIds);
dendLUT = cat(1, 0, dendLUT(:));

%% use WKW segmentation
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% get flight paths
[axonDir, axonName] = fileparts(axonFile);
axonFlightFile = fullfile(axonDir, sprintf('%s-flights.mat', axonName));

if exist(axonFlightFile, 'file')
    disp('Loading flights from cache...');
    load(axonFlightFile);
else
   [axonFlights, axonFlightsMeta] = ...
        Superagglos.getFlightPathSegIds(param, axons, 3);
    Util.saveStruct(axonFlightFile, struct( ...
        'axonFlights', axonFlights, ...
        'axonFlightsMeta', axonFlightsMeta));
end

%% only consider doubly attached flights
validFlights = axonFlightsMeta( ...
    axonFlightsMeta.numAttachments > 1, :);
validMask = ismember( ...
    axonFlights(:, {'aggloId', 'flightId'}), ...
    validFlights(:, {'aggloId', 'flightId'}), 'rows');

fprintf('\n');
fprintf('# flight paths: %d\n', numel(validMask));
fprintf('# non-dangling flight paths: %d\n', sum(validMask));

axonFlights = axonFlights(validMask, :);

%% determine nodes per flight path
[uniFlights, ~, flightNodeCount] = unique( ...
    axonFlights(:, {'aggloId', 'flightId'}), 'rows');
uniFlights.nodes = accumarray(flightNodeCount, 1);

%% accumulate evidence
axonFlights.dendId = dendLUT(1 + axonFlights.segId);
axonFlights(~axonFlights.dendId, :) = [];

[dendOverlap, ~, dendEvidence] = unique( ...
    axonFlights(:, {'aggloId', 'flightId', 'dendId'}), 'rows');
dendOverlap.evidence = accumarray(dendEvidence, 1);

% determine volume fraction of dendrite agglomerate which can be
% "explained" the the axonal flight path passing through it
dendOverlap.dendVol = accumarray( ...
    dendEvidence, axonFlights.segId, [], ...
    @(ids) sum(segSizes(unique(ids))));
dendOverlap.dendFrac = ...
    dendOverlap.dendVol ...
 ./ dendVols(dendOverlap.dendId);

% discard overlaps below evidence threshold
dendOverlap(dendOverlap.evidence < minEvidence, :) = [];

[~, rows] = ismember( ...
    dendOverlap(:, {'aggloId', 'flightId'}), ...
    uniFlights(:, {'aggloId', 'flightId'}), 'rows');
dendOverlap.flightNodes = uniFlights.nodes(rows);
dendOverlap.flightFrac = ...
    dendOverlap.evidence ...
 ./ dendOverlap.flightNodes;

% assign to axon with largest evidence
dendOverlap = sortrows(dendOverlap, 'evidence', 'descend');
[~, uniRows] = unique(dendOverlap.dendId, 'stable');
dendOverlap = dendOverlap(uniRows, :);

%%
fprintf('\n');
fprintf('# axons: %d\n', numel(axons));
fprintf('# dendrites: %d\n', numel(dends));

fprintf('\n');
fprintf('Node evidence threshold: %d\n', minEvidence);
fprintf('# dendrites with axon overlap: %d\n', size(dendOverlap, 1));

%% export examples
outputDir = '/home/amotta/Desktop/nmls';
mkdir(outputDir);

%{
% top ten
rows = 1:10;
%}

% random examples
rng(1);
rows = randperm(size(dendOverlap, 1));
rows = reshape(rows(1:10), 1, []);

%{
% good agreement between flight path and dendrite
rng(0);
rows = find( ...
    dendOverlap.flightFrac > 0.8 ...
  & dendOverlap.dendFrac > 0.8);
rows = rows(randperm(numel(rows)));
rows = reshape(rows(1:10), 1, []);
%}

for curRow = rows
    curAxonId = dendOverlap.aggloId(curRow);
    curAxon = axons(curAxonId);
    
    curDendId = dendOverlap.dendId(curRow);
    curDend = dends(curDendId);
    
    skel = skeleton();
    skel = skel.addTree( ...
        sprintf('Axon #%d', curAxonId), ...
        curAxon.nodes(:, 1:3), curAxon.edges);
    skel = skel.addTree( ...
        sprintf('Dendrite #%d', curDendId), ...
        curDend.nodes(:, 1:3), curDend.edges);
    
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel.write(fullfile(outputDir, ...
        sprintf('axon-dendrite-overlap-%d.nml', curRow)));
end