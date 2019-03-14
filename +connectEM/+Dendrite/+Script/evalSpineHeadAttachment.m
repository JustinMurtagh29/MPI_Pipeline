% Export a random subset of attached spine heads, so that we can proofread
% the quality of the final spine head attachment.
%
% Written by
%   <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
dendriteFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
outDir = '/home/amotta/Desktop';

info  = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

dendrites = load(dendriteFile);
shAgglos = dendrites.shAgglos;

dendAgglos = dendrites.dendAgglos(dendrites.indBigDends);
dendrites = dendrites.dendrites(dendrites.indBigDends);

%% Check spine head attachment
dendLUT = Agglo.buildLUT(maxSegId, dendAgglos);

shT = table;
shT.id = reshape(1:numel(shAgglos), [], 1);
shT.agglo = shAgglos(:);
shT.dendIds = cellfun( ...
    @(ids) setdiff(dendLUT(ids), 0), ...
    shT.agglo, 'UniformOutput', false);

%% Build skeletons
rng(0);

curShIds = find(cellfun(@isscalar, shT.dendIds));
curShIds = curShIds(randperm(numel(curShIds)));
curShIds = curShIds(1:25);

curDigits = ceil(log10(1 + numel(curShIds)));

for curIdx = 1:numel(curShIds)
    curShId = curShIds(curIdx);
    curShSegIds = shT.agglo{curShId};
    curDendId = shT.dendIds{curShId};
    curDend = dendrites(curDendId);
    
    curSkel = skeleton();
    curSkel = curSkel.addTree( ...
        sprintf('Dendrite %d', curDendId), ...
        curDend.nodes(:, 1:3), curDend.edges);
    curSkel = Skeleton.fromMST( ...
        segPoints(curShSegIds, :), param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf('Spine head %d', curShId);
    curSkel.colors{end} = [1, 1, 0, 1];
    
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);
    
    curOutFile = fullfile(outDir, sprintf( ...
        '%0*d_spine-head-%d.nml', curDigits, curIdx, curShId));
    curSkel.write(curOutFile);
end
