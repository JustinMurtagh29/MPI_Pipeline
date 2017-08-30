function skel = toSkel( sagglos, skel )
%TOSKEL Write superagglos to a skeleton.
% INPUT sagglos: [Nx1] struct
%           Superagglo struct array.
%       skel: (Optional) skeleton object
%           The agglos are added as trees to the provided skeleton.
%           (Default: a new skeleton is created)
% OUTPUT skel: skeleton object
%           Skeleton object containing the sagglos as trees.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('skel', 'var') || isempty(skel)
    skel = Skeleton.setParams4Dataset([], 'ex145_ROI2017');
end

for i = 1:length(sagglos)
    skel = skel.addTree(sprintf('Agglo%3d', i), sagglos(i).nodes(:,1:3), ...
        sagglos(i).edges);
end

end

