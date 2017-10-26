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

clear;
info = Util.gitInfo();

% configuration
minExplainedFrac = 0.9;
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'chiasmataSplitting/20171009T193744-kmb-on-axons-6c/outputs/20171023T115300_results.mat');
dendFile = fullfile(rootDir, 'aggloState/dendrites_03_v2.mat');
outFile = fullfile(rootDir, 'aggloState/dendrites_flight_01.mat');

%% load input data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segSizes = Seg.Global.getSegToSizeMap(param);

axons = load(axonFile);
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

dendData = load(dendFile);
dendIds = find(dendData.indBigDends);
dends = dendData.dendrites(dendIds);

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
validMask = axonFlightsMeta.numAttachments > 1;
validFlights = axonFlightsMeta(validMask, :);

fprintf('\n');
fprintf('# flight paths: %d\n', numel(validMask));
fprintf('# non-dangling flight paths: %d\n', sum(validMask));

axonFlights = axonFlights(ismember( ...
    axonFlights(:, {'aggloId', 'flightId'}), ...
    validFlights(:, {'aggloId', 'flightId'}), 'rows'), :);

%% build axon agglomerates from flight path + pickup
axonFlightAgglos = accumarray( ...
    axonFlights.aggloId, axonFlights.segId, [numel(axonIds), 1], ...
    @(ids) {setdiff(ids, 0)}, {zeros(0, 1)});

axonAgglos = Superagglos.getSegIds(axons);
axonAgglos = cellfun( ...
    @(one, two) unique(vertcat(one, two)), ...
    axonAgglos, axonFlightAgglos, 'UniformOutput', false);

axonLUT = Agglo.buildLUT(maxSegId, axonAgglos);

%% find dendrites largely explained by dendrites
dendOverlap = table;
dendOverlap.dendVol = dendVols;
dendOverlap.explainedVol = cellfun(@(ids) sum(...
    double(axonLUT(ids) > 0) .* segSizes(ids)), dendSegIds);
dendOverlap.explainedFrac = ...
    dendOverlap.explainedVol ...
 ./ dendOverlap.dendVol;

%% export examples to webKNOSSOS
%{
outputDir = '/home/amotta/Desktop/nmls';
mkdir(outputDir);

rng(0);
rows = find(...
    dendOverlap.explainedFrac < 0.9 ...
  & dendOverlap.explainedFrac > 0.5);
rows = rows(randperm(numel(rows)));
rows = reshape(rows, 1, []);

for curRow = rows(1:10)
    curDend = dends(curRow);
    curDendSegIds = dendSegIds{curRow};
    
    curAxonId = axonLUT(curDendSegIds);
    curAxonId = mode(curAxonId(curAxonId ~= 0));
    curAxon = axons(curAxonId);
    
    skel = skeleton();
    skel = skel.addTree('Axon', curAxon.nodes(:, 1:3), curAxon.edges);
    skel = skel.addTree('Dendrite', curDend.nodes(:, 1:3), curDend.edges);
    skel = Skeleton.setParams4Pipeline(skel, param);
    
    skel.write(fullfile(outputDir, ...
        sprintf('axon-dendrite-overlap-%d.nml', curRow)));
end
%}

%% evaluation
dendMask = dendOverlap.explainedFrac > minExplainedFrac;
dendOverlapIds = dendIds(dendMask);

fprintf('\n');
fprintf('# dendrite agglomerates: %d\n', numel(dendMask));
fprintf('# dendrite agglomerates explained by axon: %d\n', sum(dendMask));
fprintf('⇒ removing these\n');

%% building output
% keep all agglos except marked ones
keepMask = true(size(dendData.dendrites));
keepMask(dendOverlapIds) = false;

out = struct;
out.dendrites = dendData.dendrites(keepMask);
out.indBigDends = dendData.indBigDends(keepMask);
out.myelinDend = dendData.myelinDend(keepMask);

out.info = info;
Util.saveStruct(outFile, out);
system(sprintf('chmod a-w "%d"', outFile));