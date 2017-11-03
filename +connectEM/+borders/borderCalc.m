function borderCalc(idx)
blockSize = [32, 32, 16];

% TODO(amotta): Pass in as shared parameters (or at least paths to file)
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat', 'p');
load('/gaba/u/mberning/results/pipeline/20170217_ROI/connectomeState/connectome.mat','axons', 'dendrites');

maxSegId = Seg.Global.getMaxSegId(p);
axonLookup = Agglo.buildLUT(maxSegId, axons);
dendriteLookup = Agglo.buildLUT(maxSegId, dendrites);

% TODO(amotta):
% What about overlaps between axon and dendrite agglomerates? For now,
% axons dominate over dendrites.

% axon ids will be negative
generalLookup = dendriteLookup;
generalLookup(axonLookup ~= 0) = -axonLookup(axonLookup ~= 0);

% load segmentation and convert do double
seg = loadSegDataGlobal(p.seg, p.local(idx).bboxSmall);
blockCount = ceil(size(seg) ./ blockSize);


% apply equivalence class mapping
seg = double(seg);
segTemp = seg;
seg(seg ~= 0) = generalLookup(seg(seg ~= 0));

% set segmentation voxels that are not part of the dendrites and agglos (glia mostly) to an unused nonzero number
seg(seg == 0 & segTemp ~= 0) = Inf;
clear segTemp

% find edges and borders
% TODO(amotta): Cube borders are not properly handled yet.
[edges, ind] = connectEM.borders.codeBenedikt(seg);

corrSegId = max(seg(:)) + 1;
seg = padarray(seg, [1, 1, 1], corrSegId);
[X, Y, Z] = ind2sub(size(seg) , ind);

% correct for padding
X = X - 1;
Y = Y - 1;
Z = Z - 1;

%the +2 and -1 is for fixing the padding
blkIdx = sub2ind( ...
    blockCount, ...
    ceil(X / blockSize(1)), ...
    ceil(Y / blockSize(2)), ...
    ceil(Z / blockSize(3)));

findings = cell(prod(blockCount), 1);
areaM = cell(prod(blockCount), 1);

for curBlkIdx = 1:prod(blockCount)
    curBlkMask = (blkIdx == curBlkIdx);
    curEdges = edges(curBlkMask, :);
    
    if isempty(curEdges)
        findings{curBlkIdx} = zeros(0, 2);
        areaM{curBlkIdx} = zeros(0, 1);
        continue;
    end
    
    % count voxels per edges
   [curEdges, ~, curEdgeId] = unique(curEdges, 'rows');
    curVxCount = accumarray(curEdgeId, 1);
    
    findings{curBlkIdx} = [curEdges, curVxCount];
    
    % calculate physical area
    curVxIds = ind(curBlkMask);
    curVxIds = arrayfun( ...
        @(idx) curVxIds(curEdgeId == idx), ...
        1:size(curEdges, 1), 'UniformOutput', false);
    
    areaM{curBlkIdx} = Seg.Local.physicalBorderArea2( ...
        curVxIds, curEdges, seg, p.raw.voxelSize, true);
end

save(['/tmpscratch/kboerg/borders4/borders_' num2str(idx)],'findings','areaM');
end
