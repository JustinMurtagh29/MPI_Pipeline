function fp = getFlightPath( sagglos, mode )
%GETFLIGHTPATH Extract the flight paths from superagglos.
% INPUT sagglos: [Nx1] struct
%           Supperagglomerate struct array.
%       mode: (Optional) string
%           The mode of flight path retrieval
%           'full': (default) retrieve flight path edges and nodes
%           'nodes': retrieve only the nodes
% OUTPUT fp: struct
%           Struct in the same format as the superagglos containing only
%           the nodes and edges of flight paths. Only edges between two
%           flight path nodes are kept.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('mode', 'var') || isempty(mode)
    mode = 'full';
end

switch mode
    case 'full'
        fp = struct('nodes', [], 'edges', []);
        for i = 1:length(sagglos)
            fp(i, 1) = fpsFromSagglo(sagglos(i));
        end
    case 'nodes'
        nodes = cellfun(@(x)x(isnan(x(:,4)), 1:3), {sagglos.nodes}', ...
            'uni', 0);
        fp = cell2struct(nodes, 'nodes', 2);
    otherwise
        errro('Unknown mode %s.', mode);
end

end


function fp = fpsFromSagglo(sagglo)
idx = find(isnan(sagglo.nodes(:,4)));
fp.nodes = sagglo.nodes(idx, :);
fp.edges = sagglo.edges(all(ismember(sagglo.edges, idx), 2),:);
[~, ~, fp.edges] = unique(fp.edges);
fp.edges = reshape(fp.edges, [], 2);
end