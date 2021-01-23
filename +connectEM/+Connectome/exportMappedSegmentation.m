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
Util.clear(otherAxons, otherDends); % To reduce RAM usage
otherAgglos = Agglo.removeSegmentOverlap(otherAgglos);

allAgglos = [connAgglos(:); otherAgglos(:)];
Util.clear(connAgglos, otherAgglos);

%% Build mapping
segLUT = 1:uint64(numel(allAgglos));
segLUT = Agglo.buildLUT(maxSegId, allAgglos, segLUT);

% Build agglomerate file with mapping
curAggloFile = sprintf('%s_mapping.hdf5', runId);
curAggloFile = fullfile(outDir, curAggloFile);
curMapping = [0; segLUT(:)];

h5create( ...
    curAggloFile, '/segment_to_agglomerate', ...
    numel(curMapping), 'Datatype', class(curMapping));
h5write( ...
    curAggloFile, '/segment_to_agglomerate', curMapping);
Util.protect(curAggloFile);

% TODO(amotta): Write mapped segmentation
