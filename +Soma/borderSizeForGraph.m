function [ borderSize ] = borderSizeForGraph( borderSize, borderIdx, repC )
%BORDERSITEFORGRAPH Create a borderSize for the correspondence graph.
% INPUT borderSize: [Nx1] int/float
%           The border com list containing borders within local cubes.
%       borderIdx: [Nx1] double
%           The correspondence graph border idx array.
%       repC: (Optional) [1x1]
%           The value that is filled in for the correspondence borders.
%           (Default: NaN)
% OUTPUT borderSize: [Nx1] int/float
%           The updated border com for the correspondence graph.
% Author: Robin Hesse, based on borderComForGraph.m from Benedikt Staffler

borderIdx(isnan(borderIdx)) = size(borderSize, 1) + 1;

if ~exist('repC', 'var') || isempty(repC)
    repC = [NaN];
end
repC = repC(:)';

borderSize(end+1, :) = repC;
borderSize = borderSize(borderIdx,:);


end

