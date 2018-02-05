% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');

%% loading data
conn = load(connFile);

param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

segSizes = Seg.Global.getSegToSizeMap(param);
segCentroids = Seg.Global.getSegToCentroidMap(param);
segCentroids = segCentroids .* param.raw.voxelSize;

segCov = load(fullfile(rootDir, 'globalSegmentPCA.mat'), 'covMat');
segCov = reshape(segCov.covMat, [], 3, 3);

%% testing
axonId = find(conn.axonMeta.synCount >= 10);

rng(0);
axonId = axonId(randperm(numel(axonId), 1));
segIds = conn.axons{axonId};

segToSegDist = squareform(pdist(segCentroids(segIds, :)));
segToSegDist = segToSegDist <= 750;

nhoods = nan(sum(segToSegDist(:)), 2);
[nhoods(:, 1), nhoods(:, 2)] = find(segToSegDist);

nhoods = accumarray( ...
    nhoods(:, 1), nhoods(:, 2), ...
   [size(segToSegDist, 1), 1], @(ids) {ids});

[massesOut, comVecsOut, covMatsOut] = Agglo.mergeStatisticalMoments( ...
    segSizes(segIds), segCentroids(segIds, :), segCov(segIds, :, :), nhoods);

skel = Skeleton.fromMST( ...
    ceil(segCentroids(segIds, :) ./ param.raw.voxelSize), param.raw.voxelSize);

for curIdx = 1:numel(segIds)
    curPos = ceil(comVecsOut(curIdx, :) ./ param.raw.voxelSize);
    curCov = shiftdim(covMatsOut(curIdx, :, :));
    curCov = curCov * massesOut(curIdx);
    
    curInMat = eye(size(curCov)) * trace(curCov) - curCov;
   [~, curD] = eig(curInMat, 'vector');
   
    curAbc = sqrt(5 .* (sum(curD) / 2 - curD) ./ massesOut(curIdx));
    curAbc = ceil(repelem(curAbc(end), 3, 1) ./ param.raw.voxelSize(:))
    
    skel = skel.addTree('X', curPos + curAbc(1) .* [-1; +1] .* [1, 0, 0]);
    skel = skel.addTree('Y', curPos + curAbc(2) .* [-1; +1] .* [0, 1, 0]);
    skel = skel.addTree('Z', curPos + curAbc(3) .* [-1; +1] .* [0, 0, 1]);
end

skel = Skeleton.setParams4Pipeline(skel, param);
skel.write('/home/amotta/Desktop/blah.nml');

%{
curCov = reshape(covMatsOut, 3, 3);

%{
curWeights = segSizes(curSegIds);
curWeights = curWeights ./ sum(curWeights);

% calculate merged covariance matrix
% curCentroids = sum(curWeights .* segCentroids(curSegIds, :), 1);
curCov = sum((curWeights .* curWeights) .* segCov(curSegIds, :), 1);
curCov = reshape(curCov, 3, 3);
%}
    
curInMat = eye(size(curCov)) * trace(curCov) - curCov;
[~, curD] = eig(curInMat, 'vector');
curAbc = 2 * sqrt(5 .* (sum(curD) / 2 - curD)) ./ 1E3
%}