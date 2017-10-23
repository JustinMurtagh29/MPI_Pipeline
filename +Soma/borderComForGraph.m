function [ borderCom ] = borderComForGraph( borderCom, borderIdx, repC )
%BORDERCOMFORGRAPH Create a borderCom for the correspondence graph.
% INPUT borderCom: [Nx3] int/float
%           The border com list containing borders within local cubes.
%       borderIdx: [Nx1] double
%           The correspondence graph border idx array.
%       repC: (Optional) [1x3]
%           The value that is filled in for the correspondence borders.
%           (Default: NaN)
% OUTPUT borderCom: [Nx3] int/float
%           The updated border com for the correspondence graph.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

borderIdx(isnan(borderIdx)) = size(borderCom, 1) + 1;

if ~exist('repC', 'var') || isempty(repC)
    repC = [NaN, NaN, NaN];
elseif length(repC) == 1
    repC = repmat(repC, 1, 3);
end
repC = repC(:)';

borderCom(end+1, :) = repC;
borderCom = borderCom(borderIdx,:);


end

