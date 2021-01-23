% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

outDir = '/tmpscratch/amotta/l4/2021-01-23-mapped-segmentation';
runId = datestr(now, 30);

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

Util.log('Loading connectome');
conn = connectEM.Connectome.load(param, connFile);

Util.log('Loading axon agglomerates');
allAxons = Util.load(conn.info.param.axonFile, 'axons');
allAxons = Agglo.fromSuperAgglo(allAxons(:), true);

Util.log('Loading dendrite agglomerates');
allDends = Util.load(conn.info.param.dendriteFile, 'dendrites');
allDends = Agglo.fromSuperAgglo(allDends(:), true);

%% Building mapping
Util.log('Building mapping');
connAgglos = [conn.axons(:); conn.dendrites(:)];
connAgglos = Agglo.removeSegmentOverlap(connAgglos);

curLUT = true(maxSegId, 1);
curLUT(cell2mat(connAgglos)) = false;

otherAxons = allAxons;
otherAxons(conn.axonMeta.parentId) = [];
otherAxons = cellfun( ...
    @(segIds) segIds(curLUT(segIds)), ...
    otherAxons, 'UniformOutput', false);

otherDends = allDends;
otherDends(conn.denMeta.parentId) = [];
otherDends = cellfun( ...
    @(segIds) segIds(curLUT(segIds)), ...
    otherDends, 'UniformOutput', false);

Util.clear(curLUT);

otherAgglos = [otherAxons(:); otherDends(:)];
Util.clear(otherAxons, otherDends);
otherAgglos = Agglo.removeSegmentOverlap(otherAgglos);

allAgglos = [connAgglos(:); otherAgglos(:)];
Util.clear(connAgglos, otherAgglos);
allAgglos(cellfun(@isempty, allAgglos)) = [];

%% Build mapping
aggloIds = cellfun(@min, allAgglos);
assert(numel(aggloIds) == numel(unique(aggloIds)));
aggloIds = cast(aggloIds, 'uint64');

segLUT = 1:uint64(maxSegId);
segLUT(cell2mat(allAgglos)) = repelem( ...
    aggloIds, cellfun(@numel, allAgglos));
segLUT = [0; segLUT(:)];

% Build agglomerate file with mapping
aggloFile = fullfile(outDir, mapping.hdf5);
Util.log('Writing agglomerate file:\n%s', aggloFile);

h5create( ...
    aggloFile, '/segment_to_agglomerate', ...
    numel(segLUT), 'Datatype', class(segLUT));
h5write( ...
    aggloFile, '/segment_to_agglomerate', segLUT);
Util.protect(aggloFile);
Util.clear(segLUT);

%% Build mapped segmentation
mappedSegParam = struct;
mappedSegParam.root = fullfile(outDir, 'segmentation', '1');
mappedSegParam.backend = 'wkwrap';
mappedSegParam.dtype = 'uint32';

curTmpDir = tempname(outDir);
wkwInit('new', curTmpDir, 32, 32, mappedSegParam.dtype, 1);
wkwInit('compress', curTmpDir, mappedSegParam.root);
assert(rmdir(curTmpDir, 's'));

Util.log('Building mapped segmentation');
Seg.Global.applyMappingToSegmentation(param, aggloFile, mappedSegParam);
