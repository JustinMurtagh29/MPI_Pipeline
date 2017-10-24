function [ nodes, name, segIds ] = getSomaNodes( p, graph, meta, rp, gb, probT, sizeT, somaID )
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
% Note (BS): Why is this defined?
% mag4bbox = (mag1bbox - 1) ./4 + 1;
% mag4bbox(:,1) = ceil(mag4bbox(:,1));
% mag4bbox(:,2) = floor(mag4bbox(:,2));

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
    
% % NOTE (BS): Why separate file? Delete if possible.
% m = load('/gaba/u/mberning/results/pipeline/20170217_ROI/soma/manuNucleiMapping.mat', 'mapping');
% mapping = m.mapping;
    
%% get all the seed ids from nuclei and mapping

somaSeg = Seg.IO.loadSeg(p, bbox);

% get nucleus seg ids with at least 50% overlap with nucleus
nucleusSegIds = tabulate([0; somaSeg(nucleus)]); % add 0 just to make sure
nucleusSegIds = nucleusSegIds(2:end, 1:2);
toKeep = nucleusSegIds(:,2) > 0.5*meta.voxelCount(nucleusSegIds(:,1));
nucleusSegIds = nucleusSegIds(toKeep, 1);
clear nucleus somaSeg

% % NOTE (BS): Seems weird. Should be deleted. Although I can't fully
% % reproduce what the mapping is (e.g. for somaID 93)
% somaSeg = readKnossosRoi(p.seg.root,'2012-09-28_ex145_07x2_ROI2016_corrected_mag1', ...
%             [bbox(1,1) bbox(1,2); bbox(2,1) bbox(2,2); bbox(3,1) bbox(3,2)], 'uint32', '', 'raw');
% ind = find(nucleus);
% [x,y,z] = ind2sub(size(nucleus),ind);
% maskInd = cat(2, x, y, z);
% nucleusSegIds = unique(somaSeg(maskInd));
% nucleusSegIds = cat(1,nucleusSegIds, mapping{somaID});
   
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
nucleusAndNeighborSegIds = unique(reshape(...
    theseEdges,[2*size(theseEdges,1),1]));
center = bsxfun(@plus,4.* rp(somaID).Centroid, mag1bbox(:,1)');
% add correspondences to borderCom and borderSize
borderComWithCubeCorr = Soma.borderComForGraph...
    (gb.borderCoM, graph.borderIdx, center);
borderSizeWithCubeCorr = Soma.borderSizeForGraph...
    (gb.borderSize, graph.borderIdx, sizeT+1);
aggloSegIds = Soma.agglomerateSeedRestricted(nucleusAndNeighborSegIds,....
    center, graph.edges, meta.point, borderSizeWithCubeCorr,...
    graph.borderIdx, graph.prob, probT, sizeT,...
    borderComWithCubeCorr, somaCenterT, p.raw.voxelSize);
myEdges = graph.edges;
    
%% add all sets of seg Ids that are completely surrounded by the agglo

% get the connected components for all edges not having any agglo nodes
% and discard the largest components (= everything outside the soma)
surridx = any(ismember(myEdges, aggloSegIds),2);
theseEdges = myEdges(~surridx,:);
CCcut = Graph.findConnectedComponents(theseEdges, false, true);
maxsize = cellfun(@(C) size(C,1), CCcut);
CCcut(maxsize==max(maxsize)) = [];
aggloSegIds = sort(cat(1, aggloSegIds, cell2mat(CCcut)));

%% remove not connected components
    
% % NOTE (BS): Why is this needed?
% 
% newSegIds = unique(aggloSegIds);
% surridx = any(xor(ismember(graph.edges, newSegIds),2),2);
% theseEdges = graph.edges(~surridx,:);
% CCcut = Graph.findConnectedComponents(theseEdges, false, true);
% 
% if size(CCcut,1) > 1
%     [maxsize] = cellfun(@(C) size(C,1), CCcut);
%     %just to test if connected component is somatic
%     if max(maxsize) < 10000 
%         CCcut(maxsize~=max(maxsize)) = [];
%     end
% end
% aggloSegIds = CCcut{1};

%% return
aggloSegIds=unique(aggloSegIds);
segIds = aggloSegIds;
nodes = meta.point(:,aggloSegIds)';
name = sprintf('soma%d_%.2f_%d_Agglo', somaID, probT, sizeT);
    
end



