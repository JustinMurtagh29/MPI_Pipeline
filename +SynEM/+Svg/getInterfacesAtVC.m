function [interfaces, edges, borders] = getInterfacesAtVC(seg, vcSeg, r, voxelSize)
%GETINTERFACESATVC Calculate segment interfaces  close to vesicle clouds.
% INPUT seg: 3d int array
%           Volume segmentation.
%       vcSeg: 3d logical array
%           Vesicle segmentation for the same region as seg.
%       r: (Optional) double
%           Size of the spherical structuring element used for dilation of
%           vcSeg in units of voxel size.
%           (Default: No dilation of vcSeg)
%       voxelSize: (Optional) [1x3] double
%           Voxel size of seg/vcSeg.
%           (Default: [1 1 1])
% OUTPUT edges: [Nx2] int
%           Integer ids for VC edges between segments.
%           (see also SynEM.Svg.findEdgesAndBorders)
%        borders: [Nx1] struct
%           VC border struct array.
%           (see also SynEM.Svg.findEdgesAndBorders)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if exist('r','var') && ~isempty(r)
    if ~exist('voxelSize','var')
        voxelSize = [];
    end
    h = Util.ballMask(r, voxelSize);
    vcSeg = imdilate(vcSeg, h);
end

[edges, borders] = SynEM.Svg.findEdgesAndBorders(seg, ~vcSeg);
[interfaces, intIdx] = SynEM.Svg.calculateInterfaces(seg, ...
    edges, borders, 150, voxelSize, [40 80 160]);
edges = edges(intIdx,:);
borders = borders(intIdx,:);
end
