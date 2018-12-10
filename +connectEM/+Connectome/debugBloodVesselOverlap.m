% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

outDir = '/home/amotta/Desktop';

vesselScoreThresh = 0.5;
vesselFracThresh = 0.05;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

vessel = fullfile(rootDir, 'heuristicResult.mat');
vessel = load(vessel);
vessel.segIds = double(vessel.segIds);
vessel = vessel.segIds(vessel.vesselScore > vesselScoreThresh);

maxSegId = Seg.Global.getMaxSegId(param);
vesselLUT = logical(Agglo.buildLUT(maxSegId, {vessel}));

segPoints = Seg.Global.getSegToPointMap(param);
segSizes = Seg.Global.getSegToSizeMap(param);

%% Calculate blood vessel fraction per axon
axonBvFrac = cellfun(@(ids) ...
    sum(vesselLUT(ids) .* segSizes(ids)) ...
    / sum(segSizes(ids)), conn.axons);
axonIds = find(axonBvFrac > vesselFracThresh);

%% Report fraction of volume and path length
% TODO(amotta): Verify that this path length measure is proportional to the
% one used for the numbers reported in the manuscript.
axonVols = cellfun(@(ids) sum(segSizes(ids)), conn.axons);
axonLens = conn.axonMeta.pathLen;

bvFrac = axonVols / sum(axonVols);
bvFrac = sum(axonBvFrac .* bvFrac) %#ok
bvAxonVolFrac = sum(axonVols(axonIds)) / sum(axonVols) %#ok
bvAxonPathLenFrac = sum(axonLens(axonIds)) / sum(axonLens) %#ok

%% Export to webKnossos
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(axonIds)
    curAxonId = axonIds(curIdx);
    curAxonSegIds = conn.axons{curAxonId};
    
    curSkel = skel;
    curSkel = Skeleton.fromMST( ...
        segPoints(curAxonSegIds, :), param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf('Axon %d', curAxonId);
    curSkel.colors{end} = [0, 0, 1, 1];
    
    curSkelName = sprintf( ...
        '%0*d_axon-%d.nml', ...
        ceil(log10(1 + numel(axonIds))), curIdx, curAxonId);
    curSkel.write(fullfile(outDir, curSkelName));
end
