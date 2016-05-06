function [ hasOv ] = isInTrainingDataBBox( stacks, bbox )
%ISINTRAININGDATABBOX Check whether the a specified bbox overlaps with the
%cortex training data regions.
% INPUT stacks: Cortex training data parameter struct.
%       bbox: [3x2] array of integer specifying a rectangular bounding box
%           in x, y and z dimension.
% OUTPUT hasOv: [Nx1] of logical containing true if the corresponding
%           stack and bbox have any overlap.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

hasOv = false(length(stacks),1);

for i = 1:length(stacks)
    stackBBox = stacks(i).bboxRaw;
    hasOv(i) = all(bsxfun(@le,bbox(:,1),stackBBox(:,2))) & ...
               all(bsxfun(@ge,bbox(:,2),stackBBox(:,1)));
end

end

