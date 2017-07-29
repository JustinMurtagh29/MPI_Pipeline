function [interfaces, edges, borders] = getInterfacesAtVC(seg, vcSeg, ...
    r, voxelSize, r_merge)
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
%       r_merge: (Optional) double
%           Size of the structuring element used for erosion of the vcSeg.
%           All segments that overlap with one eroded connected components
%           will be merged in the segmentation. Set 0 to do the segment
%           merging without prior erosion.
%           (Default: no segment merging).
% OUTPUT edges: [Nx2] int
%           Integer ids for VC edges between segments.
%           (see also SynEM.Svg.findEdgesAndBorders)
%        borders: [Nx1] struct
%           VC border struct array.
%           (see also SynEM.Svg.findEdgesAndBorders)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if exist('r_merge', 'var') && ~isempty(r_merge)
    if r_merge > 0
        vcOv = imerode(vcSeg, Util.ballMask(30, voxelSize));
    end
    stats = regionprops(vcOv, 'PixelIdxList');
    for i = 1:length(stats)
        toComb = setdiff(seg(stats(i).PixelIdxList), 0);
        if ~isempty(toComb)
            seg(ismember(seg,toComb)) = toComb(1);
        end
    end
    seg = Seg.Local.closeGaps(seg);
end

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
