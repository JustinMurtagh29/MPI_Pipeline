function [ interfaceSurfaceList, centroidList, neighborIDs, gpProb ] = loadInterfacesFromSeg( cube, areaThreshold, idConversion )
%LOADINTERFACESFROMSEG Load borders, centroids and edges.
% INPUT cube: Parameter file of a seg cube (p.local(i)).
%       areaThreshold: Consider only interfaces with > areaThreshold
%                      voxels.
%       idConversion: String specifying ids conversion of edges. Currently
%                     the options are 'none','globalToLocal'.
%                     (Standard:none)
% OUTPUT interfaceSurfaceList: List of interfaces linear voxels in cube
%                     bboxSmall
%        centroidList: Cell array with centroids of corresponding
%                     interface.
%        neighborIDs: Matrix containing the two adjacent segments of each
%                     interface.
%        gpProb: (Gaussian process) probability of the two edges to belong
%                together.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('idConversion','var') || isempty(idConversion)
    idConversion = 'none';
end

%load borders, edges and centroids
m = matfile(cube.borderFile);
borders = m.borders;
m = matfile(cube.edgeFile);
edges = m.edges;
borders = squeeze(struct2cell(borders))';
keepIndices = cellfun(@(x) x > areaThreshold, borders(:,2));
interfaceSurfaceList = borders(keepIndices,1);
interfaceSurfaceList = cellfun(@(x)double(x'),interfaceSurfaceList,'UniformOutput',false);
centroidList = borders(keepIndices,3);
neighborIDs = edges(keepIndices,:);
if nargout == 4
    m = matfile(cube.probFile);
    prob = m.prob;
    gpProb = prob(keepIndices);
end
if ~strcmp(idConversion,'none')
    neighborIDs = Seg.Local.localGlobalIDConversion(idConversion, cube, neighborIDs);
end

end

