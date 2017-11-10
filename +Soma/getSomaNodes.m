function [ nodes, name, segIds ] = getSomaNodes( p, graph, meta, rp, ...
    gb, probT, sizeT, somaID )
%GETSOMANODES Agglomerate soma using Christians nuclei and Manuels list
% INPUT     p, graph, meta, rp, gb: connectomics stuff 
%               default data that is used.
%           probT: int
%               Only edges with prob > probThreshold will be taken.
%           sizeT: int
%               Merge probability for the corresponding edges.
%           somaID: int
%               ID of the soma (corresponding to Christians nuclei)
% OUTPUT    nodes: [Nx3] int
%               meta.point of ids that are in the agglo.
%           name: str
%               The edges that connect the ids.
%           segIds: [Nx1]
%               seg ids that are in the agglo.
%           The edges that connect the ids.
% Author: Robin Hesse
% Modified by: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% mag stuff to get wk coords
mag1bbox = [128, 128, 128, 5573, 8508, 3413]; 
mag1bbox = Util.convertWebknossosToMatlabBbox(mag1bbox);
mag1bbox = mag1bbox + [-25 25; -25 25; -10 10];

margin = [15, 15, 10]; 
bboxM4Cropped = [round(rp(somaID).BoundingBox(1:3)) - margin; ...
    round(rp(somaID).BoundingBox(1:3) + rp(somaID).BoundingBox(4:6)) ...
    + margin]';
bbox = bsxfun(@plus,bsxfun(@plus,4.*bboxM4Cropped,[-1, 2]),  mag1bbox(:,1));


%% load nuclei und segId-Nuclei-Mapping
    
% Note (BS): Load directly from nuclei segmentation?
m = load(strcat(['/gaba/u/mberning/results/pipeline/20170217_ROI/' ...
    'soma/Nuclei/Nucleus'], int2str(somaID), '.mat'));
nucleus = m.nucleus;
if isfield(m, 'bbox') % load bbox from file if exists
    bbox = m.bbox;
else
    % fix of bbox from nuclear_pores/+MouseROI2016/ProcessSingleNucleus2.m
    raw = readKnossosRoi(p.raw.root, p.raw.prefix, bbox);
    if min(min(min(raw))) == 0
        boxsize = size(raw)';
        mask = raw==0;
        maskrp = regionprops(permute(~mask, [2 1 3]));
        maskbox = [ceil(maskrp.BoundingBox(1:3))', ...
            floor(maskrp.BoundingBox(1:3) + maskrp.BoundingBox(4:6))'];
        bbox(:,1) = bsxfun(@plus, ...
            bsxfun(@plus, bbox(:,1), maskbox(:,1)), -1);
        bbox(:,2) = bsxfun(@plus, bbox(:,2), -(boxsize - maskbox(:,2)));
    end
    clear raw mask maskrp maskbox boxsize
end


%% get all the seed ids from nuclei and mapping

somaSeg = Seg.IO.loadSeg(p, bbox);

% this caused some bugs so assert it
assert(all(size(somaSeg) == size(nucleus)));

% get nucleus seg ids with at least 50% overlap with nucleus
nucleusSegIds = tabulate([0; somaSeg(nucleus)]); % add 0 just to make sure
nucleusSegIds = nucleusSegIds(2:end, 1:2);
toKeep = nucleusSegIds(:,2) > 0.5*meta.voxelCount(nucleusSegIds(:,1));
nucleusSegIds = nucleusSegIds(toKeep, 1);
clear nucleus somaSeg
   

%% if no seed was found fill in dummy stuff, so parallel function doesn't break
if size(nucleusSegIds,1) == 0
    segIds = nan;
    center = bsxfun(@plus,4.* rp(somaID).Centroid, mag1bbox(:,1)');
    nodes = center;
    name = sprintf('soma%d_%.2f_%d_Agglo', somaID, probT, sizeT);
    return;
end

    
%% get all neighbors of the nuclei seg Ids as seeds

somaCenterT = 8000; % maximal distance of agglos to seed point
idx = any(ismember(graph.edges, nucleusSegIds),2);
theseEdges = graph.edges(idx,:);
nucleusAndNeighborSegIds = unique(theseEdges(:));
center = bsxfun(@plus,4.* rp(somaID).Centroid, mag1bbox(:,1)');
% add correspondences to borderCom and borderSize
borderComWithCubeCorr = Soma.borderComForGraph...
    (gb.borderCoM, graph.borderIdx, center);
borderSizeWithCubeCorr = Soma.borderSizeForGraph...
    (gb.borderSize, graph.borderIdx, sizeT+1);
aggloSegIds = Soma.agglomerateSeedRestricted(nucleusAndNeighborSegIds,....
    center, graph.edges, meta.point, borderSizeWithCubeCorr,...
    graph.prob, probT, sizeT,...
    borderComWithCubeCorr, somaCenterT, p.raw.voxelSize);


%% add all sets of seg Ids that are completely surrounded by the agglo

% get the connected components for all edges not having any agglo nodes
% and discard the largest components (= everything outside the soma)
surridx = any(ismember(graph.edges, aggloSegIds),2);
CCcut = Graph.findConnectedComponents(graph.edges(~surridx,:), false, false);
maxsize = cellfun(@(C) size(C,1), CCcut);
CCcut(maxsize==max(maxsize)) = [];
aggloSegIds = unique(cat(1, aggloSegIds, cell2mat(CCcut))); % just to make sure


%% return

segIds = aggloSegIds;
nodes = meta.point(:,aggloSegIds)';
name = sprintf('soma%d_%.2f_%d_Agglo', somaID, probT, sizeT);
    

end



