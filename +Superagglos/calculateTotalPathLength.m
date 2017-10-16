function TPL = calculateTotalPathLength(agglos,scale)
% returns for each superagglo the total path length (unit dependent on
% scale input argument)
% Caution: as this function calculates the interdistance by using the superagglo
% coordinates which basically are segment center points, this is only an 
% approximation and the more inaccurate the smaller the superagglo (with
% single segments having TPL of 0 nm)
%
% INPUT
% agglos        agglos in the superagglo format
% scale         1x3 vector with the voxel sizes in case the superagglo
%               coordinates are pixel bases (what they normally are)
% OUTPUT
% TPL           total path length for each superagglo with units the same
%               as "scale"

if ~exist('scale','var') || isempty(scale)
    scale = [1 1 1];
    warning('No scale was given as input. Assuming voxel size of 1 nm')
end


TPL = arrayfun(@(x) sum(arrayfun(@(y) sqrt(sum(diff(bsxfun(@times,x.nodes(x.edges(y,:)',1:3),scale),1,1).^2)),1:size(x.edges,1))), agglos);