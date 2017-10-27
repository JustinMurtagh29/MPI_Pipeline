function bboxes = nucleiBbox( rp, p )
%NUCLEIBBOX Get the bbox of the nuclei.
% INPUT rp: [Nx1] struct
%           Nuclei metadata struct. Needs to contain the field BoundingBox.
%       p: struct
%           Segmentation parameter struct for ex145_07x2_ROI2017
% OUTPUT bboxes: [Nx1] cell
%           Bounding boxes for the respective nuclei.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

mag1bbox = [128, 128, 128, 5573, 8508, 3413]; 
mag1bbox = Util.convertWebknossosToMatlabBbox(mag1bbox);
mag1bbox = mag1bbox + [-25 25; -25 25; -10 10];
margin = [15, 15, 10]; 

bboxes = cell(length(rp), 1);
for i = 1:length(rp)
    
    bboxM4Cropped = [round(rp(somaID).BoundingBox(1:3)) - margin; ...
        round(rp(i).BoundingBox(1:3) + rp(i).BoundingBox(4:6)) ...
        + margin]';
    bbox = bsxfun(@plus,bsxfun(@plus,4.*bboxM4Cropped,[-1, 2]), ...
        mag1bbox(:,1));
    
    % code from nuclear_pores/+MouseROI2016/ProcessSingleNucleus2.m
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
    
    bboxes{i} = bbox;
end

end

