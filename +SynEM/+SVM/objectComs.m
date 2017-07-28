function [ com ] = objectComs( seg, bbox )
%OBJECTCOMS Get the center of mass for each object in the svm segmentation.
% INPUT seg: 3d int
%           Segmentation with ids for SVM objects ('1' synapses, '2'
%           vesicle clouds, '3' mito)
%       bbox: (Optional) [3x2] or [3x1] int
%           Bounding box for seg or global coordinate of top left
%           coordinate in seg (i.e. bbox(:,1))
%           (Default: output centroids are not translated to global bbox
%           coordinates)
% OUTPUT com: struct
%           Struct containing the coms for the SVM objects.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if exist('bbox', 'var')
    bbox = bbox(1:3);
    bbox = bbox(:)';
end

cen = cell(3, 1);
for i = 1:3
    stats = regionprops(seg == i, 'Centroid');
    cen{i} = reshape([stats.Centroid], 3, [])';
    cen{i} = round(cen{i}(:, [2 1 3])); % regionprops interchanges x and y
    if exist('bbox', 'var')
        cen{i} = bsxfun(@plus, cen{i}, bbox - 1);
    end
    
end

com.synapses = cen{1};
com.vesicle = cen{2};
com.mito = cen{3};


end

