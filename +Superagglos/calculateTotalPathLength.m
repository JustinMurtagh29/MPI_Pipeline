function TPL = calculateTotalPathLength(agglos,scale)
% returns for each superagglo the total path length in nm
% Caution: as this function calculates the interdistance by using the superagglo
% coordinates which basically are segment center points, this is only an 
% approximation and the more inaccurate the smaller the superagglo (with
% single segments having TPL of 0 nm)


if ~exist('scale','var') || isempty(scale)
    scale = [1 1 1];
    warning('No scale was given as input. Assuming voxel size of 1 nm')
end


arrayfun(@(x) sum(arrayfun(@(y) sqrt(sum(diff(bsxfun(@times,x.nodes(x.edges(y,:)',1:3),scale),1,1).^2)),1:size(x.edges,1))), agglos)