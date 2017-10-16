function l = mstLength( sagglos, scale )
%MSTLENGTH Calculate the length of an agglo via its minimal spanning tree.
% INPUT agglos: struct
%           Struct array of superagglos.
%       scale: (Optional) [1x3] double
%           Voxel size scale.
%           (Default: [1, 1, 1])
% OUTPUT l: [Nx1] double
%           The MST lenght of the agglos.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('scale', 'var') || isempty(scale)
    scale = false;
end

l = zeros(length(sagglos), 1);
for i = 1:length(sagglos)
    
    % single node agglos
    if size(sagglos(i).nodes, 1) < 2; continue; end;
    
    % MST
    if scale
        distMat = sparse(squareform(pdist( ...
            bsxfun(@times, sagglos(i).nodes(:,1:3), scale(:)'))));
    else
        distMat = sparse(squareform(pdist(sagglos(i).nodes(:,1:3))));
    end
    adjMat = graphminspantree(distMat);
    l(i) = full(sum(adjMat(:)));
end

end
