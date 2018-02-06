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
    ceil( ...
        segCentroids(segIds, :) ...
        ./ param.raw.voxelSize), ...
	param.raw.voxelSize);

for curIdx = 1:numel(segIds)
    curMass = massesOut(curIdx);
    curPos = comVecsOut(curIdx, :);
    curCov = shiftdim(covMatsOut(curIdx, :, :));
    
    % convert covariance matrix to tensor of inertia
    curCov = curCov * curMass;
    curInMat = eye(size(curCov)) * trace(curCov) - curCov;
   
    % extract half-axes of ellipsoid from tensor of inertia
   [~, curD] = eig(curInMat, 'vector');
    curAbc = sqrt(5 .* (sum(curD) / 2 - curD) ./ curMass);
    curAbc = repelem(sqrt(prod(curAbc(2:3))), 3, 1);
    curAbc = ceil(curAbc ./ param.raw.voxelSize(:));
    
    curPos = ceil(curPos ./ param.raw.voxelSize);
    skel = skel.addTree('X', curPos + curAbc(1) .* [-1; +1] .* [1, 0, 0]);
    skel = skel.addTree('Y', curPos + curAbc(2) .* [-1; +1] .* [0, 1, 0]);
    skel = skel.addTree('Z', curPos + curAbc(3) .* [-1; +1] .* [0, 0, 1]);
end

skel = Skeleton.setParams4Pipeline(skel, param);

skelName = sprintf('axon-%-diameter.nml', axonId);
skel.write(fullfile('/home/amotta/Desktop', skelName));