function nodes = getEndingsDendrite(segmentMeta, agglo)
if iscell(agglo)
    if numel(agglo)>1
        error('agglo has to be single agglo')
    end
    agglo = agglo{1};
end
% get smoothed BW data (and offset)
[voxelRepresentation2, min_coord] = connectEM.findQueries(agglo, segmentMeta, struct('debug', false));

%create skeleton from BW data
if ~exist('Skeleton3D','file')
addpath('/gaba/u/kboerg/code/skeleton3d-matlab')
end
if ~exist('Skel2Graph3D','file')
    addpath('/gaba/u/kboerg/code/skel2graph3d-matlab/')
end
skel = Skeleton3D(voxelRepresentation2>0);
if  any(skel(:))
    node2 = connectEM.querySkeleton(skel);
    nodes = cell2mat(arrayfun(@(x) round([x.comx * 8, x.comy*8, x.comz*8/2.313]+min_coord),node2(arrayfun(@(x) length(x.links)==1,node2)),'uni',0));  % move again into dataset coords
else
    nodes = [];
end
end
