load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
connectome = load('/gaba/u/mberning/results/pipeline/20170217_ROI/connectomeState/connectome.mat');

% TODO(amotta): Pass in as argument
blockSize = [32, 32, 16];

cubeCount = size(p.local);
boxSize = 1 + diff(p.bbox, 1, 2)';
tileSize = reshape(p.tileSize, 1, []);

blockCount = ceil(boxSize ./ blockSize);
blocksPerCube = ceil(tileSize ./ blockSize);

globalFindings = cell(blockCount);
globalAreaM = cell(blockCount);

for curCubeIdx = 1:prod(cubeCount)
    disp(curCubeIdx);
    
    curCubeSize = diff(p.local(curCubeIdx).bboxSmall, 1, 2);
    curCubeSize = 1 + reshape(curCubeSize, 1, []);
    
    % load pre-computed per-cube surface data
    curResults = load(['/tmpscratch/kboerg/borders3/borders_' num2str(curCubeIdx)]);
    
    curFindings = curResults.findings;
    curAreaM = curResults.areaM;
    
    off = nan(1, 3);
   [off(1), off(2), off(3)] = ind2sub(cubeCount, curCubeIdx);
    off = (off - 1) .* blocksPerCube + 1;
    
    globalFindings( ...
        off(1):(off(1) + size(curFindings, 1) - 1), ...
        off(2):(off(2) + size(curFindings, 2) - 1), ...
        off(3):(off(3) + size(curFindings, 3) - 1)) = curFindings;
    globalAreaM( ...
        off(1):(off(1) + size(curAreaM, 1) - 1), ...
        off(2):(off(2) + size(curAreaM, 2) - 1), ...
        off(3):(off(3) + size(curAreaM, 3) - 1)) = curAreaM;
end

axonId = -28;
globalAxon = zeros(blockCount);
targetAll = zeros(blockCount);
targetSmooth = zeros(blockCount);
targetAD = zeros(blockCount);

idsSmooth = double(connectome.idxSD);
idsAD = double(connectome.idxAD);

for curX = 1:blockCount(1)
    disp(curX)
    for curY = 1:blockCount(2)
        for curZ = 1:blockCount(3)
            curFindings = globalFindings{curX, curY, curZ};
            curAreas = globalAreaM{curX, curY, curZ};
            curEdges = curFindings(:, 1:2);
            
            globalAxon(curX, curY, curZ) = ...
                sum(curFindings(any(ismember(curEdges, axonId), 2), 3));
            targetAll(curX, curY, curZ) =  ...
                sum(sum(curEdges > 0 & ~isinf(curEdges), 2) .* curAreas);
            
            % NOTE(amotta): If we can guarantee that target classes are
            % pairwise disjoint (which they should be), the following can
            % be achieved with a single call to `ismember`, which leads to
            % a significant speed-up.
            targetSmooth(curX, curY, curZ) = ...
                sum(sum(ismember(curEdges, idsSmooth), 2) .* curAreas);
            targetAD(curX, curY, curZ) = ...
                sum(sum(ismember(curEdges, idsAD), 2) .* curAreas);
        end
    end
end

% no imresize3 on cluster MATLAB
assert(p.raw.voxelSize(1) == p.raw.voxelSize(2));
sG = size(globalAxon);
sGS = round(sG .* [1,1, p.raw.voxelSize(3) / (2 * p.raw.voxelSize(2))]);
globalAxonStretch = zeros(sGS);

for curX = 1 : sG(1)
    globalAxonStretch(curX,:,:) = permute(imresize(squeeze(globalAxon(curX,:, :)), sGS(2:3)),[3, 1, 2]);
end

distanceT_Stretch = bwdist(globalAxonStretch>0);
distanceT = zeros(sG);
for curX = 1 : sGS(1)
    distanceT(curX, :, :) = permute(imresize(squeeze(distanceT_Stretch(curX,:,:)), sG(2:3)),[3, 1, 2]);
end

for i_dist = 0 : 60
    selector = @(x)sum(sum(sum(x(distanceT>=i_dist&distanceT<i_dist+1))));
    graph_AD(i_dist+1) = selector(targetAD)/selector(targetAll);
    graph_Smooth(i_dist+1) = selector(targetSmooth)/selector(targetAll);
end    
    