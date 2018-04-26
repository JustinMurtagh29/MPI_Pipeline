% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear all; %#ok

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_19_a.mat');

nmlStart = 401;
nmlCount = 400;
nmlFileCount = 50;
nmlDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

axons = load(axonFile);
axonIds = find(axons.indBigAxons);
axons = Agglo.fromSuperAgglo(axons.axons(axonIds));

%% Sort axons by size
[~, sortIds] = sort(cellfun( ...
    @numel, axons), 'descend');
% Ignore percolator
sortIds(1) = [];

axons = axons(sortIds);
axonIds = axonIds(sortIds);

%% Export
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

numDigits = ceil(log10(nmlStart + nmlCount));
fileIdStart = ceil(nmlStart / nmlFileCount);
fileIdEnd = ceil((nmlStart + nmlCount - 1) / nmlFileCount);

for curId = 1:numFiles
    curOff = 1 + (curId - 1) * nmlFileCount;
    curIdStart = max(curOff, nmlStart);
    curIdsEnd = min(curOff + nmlFileCount, nmlStart + nmlCount) - 1;
    curIds = curIdStart:curIdsEnd;
    
    curAxons = axons(curIds(:));
    curAxonIds = axonIds(curIds(:));
    
    curSkel = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        curAxons, 'UniformOutput', false);
    curSkel = Skeleton.fromMST( ...
        curSkel, param.raw.voxelSize, skel);
    curSkel.names = arrayfun( ...
        @(idx, id) sprintf('%0*d. Axon %d', numDigits, idx, id), ...
        curIds(:), curAxonIds, 'UniformOutput', false);
    
    curSkel.write(fullfile(nmlDir, sprintf( ...
        'largest-axons_%d-to-%d.nml', curIds(1), curIds(end))));
end
