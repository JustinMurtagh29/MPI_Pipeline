function agglo = deleteNodes(agglo, nodeIdx, keepMode)
%DELETENODES Delete the specified nodes from superagglos.
%
% INPUT agglos: [Nx1] struct
%           The superagglos.
%       nodeIdx: [Nx1] cell
%           Linear or logical indices of the nodes that are deleted from
%           the corresponding agglomerate.
%       keepMode: (Optional) logical
%           Keep only the specified nodes instead of deleting them.
%           (Default: false)
%
% OUPUT agglo: [Nx1] struct
%           The resulting superagglomerates.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('keepMode', 'var') || isempty(keepMode)
    keepMode = false;
end

assert(length(agglo) == length(nodeIdx));
agglo = Superagglos.origFields(agglo);

for i = 1:length(agglo)
    
    idx = nodeIdx{i};
    if islogical(idx)
        lidx = idx(:);
        idx = find(lidx);
    else
        lidx = false(size(agglo.nodes, 1), 1);
        lidx(idx) = true;
    end
    
    if keepMode
        lidx = ~lidx;
        idx = find(lidx);
    end

    agglo(i).nodes(lidx, :) = [];
    agglo(i).edges(any(ismember(agglo(i).edges, idx), 2), :) = [];
    decr = cumsum(lidx);
    agglo(i).edges = agglo(i).edges - decr(agglo(i).edges);
end

end

