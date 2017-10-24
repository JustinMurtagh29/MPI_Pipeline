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
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'chiasmataSplitting/20171009T193744-kmb-on-axons-6c/outputs/20171023T115300_results.mat');
dendFile = fullfile(rootDir, 'aggloState/dendrites_03_v2_splitmerged.mat');

%% load input data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

axons = load(axonFile);
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

dends = load(dendFile);
dendIds = find(dends.indBigDends);
dends = dends.dendrites(dendIds);

%% use WKW segmentation
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% get flight paths
[axonDir, axonName] = fileparts(axonFile);
axonFlightFile = fullfile(axonDir, sprintf('%s-flights.mat', axonName));

if exist(axonFlightFile, 'file')
    disp('Loading flights from cache...');
    load(axonFlightFile, 'axonFlights');
else
    axonFlights = ...
        Superagglos.getFlightPathSegIds(param, axons, 3);
    Util.save(axonFlightFile, axonFlights);
end

%% check overlap with dendrites
maxSegId = Seg.Global.getMaxSegId(param);
dendSegIds = Superagglos.getSegIds(dends);

dendLUT = Agglo.buildLUT(maxSegId, dendSegIds);
dendLUT = cat(1, 0, dendLUT(:));

%% accumulate evidence
axonFlights.dendId = dendLUT(1 + axonFlights.segId);

[axonDendOverlap, ~, axonDendCount] = unique( ...
    axonFlights(:, {'axonId', 'flightId', 'dendId'}), 'rows');
axonDendOverlap.evidence = accumarray(axonDendCount, 1);
axonDendOverlap(~axonDendOverlap.dendId, :) = [];