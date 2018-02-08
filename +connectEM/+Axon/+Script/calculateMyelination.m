% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');

myelinThresh = 0.5;

%% loading data
conn = load(connFile);

% make sure that axons are pairwise non-overlapping
axonSegIds = cat(1, conn.axons{:});
assert(numel(axonSegIds) == numel(unique(axonSegIds)));
clear axonSegIds;

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

segCentroids = Seg.Global.getSegToCentroidMap(param);
segCentroids = segCentroids .* param.raw.voxelSize;

graph = fullfile(rootDir, 'graph.mat');
graph = load(graph, 'edges', 'borderIdx');

borders = fullfile(rootDir, 'globalBorder.mat');
borders = load(borders, 'borderArea2');

heuristics = fullfile(rootDir, 'heuristicResult.mat');
heuristics = load(heuristics, 'segIds', 'myelinScore');

%% build look-up tables
myelinLUT = false(maxSegId, 1);
myelinLUT(heuristics.segIds) = ...
    heuristics.myelinScore > myelinThresh;

axonLUT = Agglo.buildLUT(maxSegId, conn.axons);

%% build graph table
graph = struct2table(graph);

% remove correspondences (i.e., non-borders)
graph(isnan(graph.borderIdx), :) = [];

% remove intra-agglomerate borders
graph.axons = axonLUT(graph.edges);
graph(graph.axons(:, 1) ~= graph.axons(:, 2), :) = [];

%% per-segment
seg = table;
seg.id = graph.edges(:);
seg.borderArea = repmat(borders.borderArea2(graph.borderIdx), 2, 1);
seg.isMyelin = repmat(any(myelinLUT(graph.edges), 2), 2, 1);

seg.myelinArea = seg.isMyelin .* seg.borderArea;
seg.isMyelin = [];

areaLUT = table;
areaLUT.border = accumarray(seg.id, seg.borderArea, [maxSegId, 1]);
areaLUT.myelin = accumarray(seg.id, seg.myelinArea, [maxSegId, 1]);
assert(all(areaLUT.myelin <= areaLUT.border));

%% along axon
% NOTE(amotta): We've now calculated the total border area and myelin
% contact area for each segment. These values could now be pooled over an
% agglomerate.
%
% In reality, however, the degree of myelination changes along an axon. So,
% let's do a "marching sphere" approach (similar to what was done for the
% diameter calculation).
nhoodDistThresh = 1500;

myelination = cell(size(conn.axons));
for curAxonIdx = 1:numel(conn.axons)
    curAxonIdx
    curSegIds = conn.axons{curAxonIdx};
    curCentroids = segCentroids(curSegIds, :);
    
    % build minimal spanning tree
    curDists = sparse(squareform(pdist(curCentroids)));
    
    %{
    % NOTE(amotta): Uncomment this block to build segment neighbourhoods
    % based on the distance along the minimum spanning tree of the
    % agglomerate. For branchpoints and parallel running processes this
    % might yield slightly better results (albeit at higher computational
    % cost).
    curDists = graphminspantree(curDists);
    curDists = graphallshortestpaths(curDists, 'Directed', false);
    %}
    
    curMyelination = nan(size(curSegIds));
    for curSegIdx = 1:numel(curSegIds)
        curNhoodSegIds = curDists(:, curSegIdx) < nhoodDistThresh;
        curNhoodSegIds = curSegIds(curNhoodSegIds);
        
        curMyelination(curSegIdx) = ...
            sum(areaLUT.myelin(curNhoodSegIds)) ...
          / sum(areaLUT.border(curNhoodSegIds));
    end
    
    myelination{curAxonIdx} = curMyelination;
end