function [ synapses, boutons ] = clusterBoutonSynapsesLocalPostsynAgglo( ...
    boutons, edges, contProb, probT, dist)
%CLUSTERBOUTONSYNAPSESLOCALPOSTSYNAGGLO Cluster the synapses of boutons
%based on a local agglomeration of the postsynaptic segments.
% INPUT boutons: [Nx3] table
%           Table of boutons agglomerations containing the columns
%           postsynId and edgeIdx.
%       edges: [Nx2] int
%           The global edge list.
%       contProb: [Nx1] double
%           The continuity probability for the corresponding edges.
%       probT: double
%           Lower threshold on the continuity probability for local
%           postsynaptic agglomeration.
%       dist: int
%           Maximal graph distance from postsynaptic segments for
%           postsynaptic agglomeration.
% OUTPUT synapses: table
%           Table containing the synapse information for the corresonding
%           bouton consisting of number of synapses (n_syn), all edge
%           indices that form each synapse (edgeIdx) and the all
%           postsynaptic ids (targetIdx).
%        boutons: table
%           The input table with the synapses information concatenated as
%           additional rows.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% get edges above threshold
edgesT = edges(contProb > probT, :);

% get the neighbor list
nIds = Graph.edges2Neighbors(edgesT);

noB = size(boutons, 1);
synapses.syn_n = zeros(noB, 1);
synapses.syn_edgeIdx = cell(noB, 1);
synapses.syn_postsynId = cell(noB, 1);
fprintf('Clustering synaptic interfaces:  0%%');
for i = 1:noB
    ps = boutons.postsynId{i};
    [ups, ~, ic] = unique(ps);
    edgeIdx = boutons.edgeIdx{i};
    if length(ups) > 1
        [~, uT] = L4.Agglo.distRestrictedAgglo(ups, nIds, [], [], dist, true);
    else
        uT = 1;
    end
    T = uT(ic);
    synapses.syn_n(i) = max(T);
    synapses.syn_edgeIdx{i} = accumarray(T, edgeIdx, [], @(x){x});
    synapses.syn_postsynId{i} = accumarray(T, ps, [], @(x){x});
    
    if floor(i/noB*100) < 10
        fprintf('\b\b%d%%',floor(i/noB*100));
    else
        fprintf('\b\b\b%d%%',floor(i/noB*100));
    end
    
end
fprintf('\n');

synapses.syn_presynId = getPresynId([synapses.syn_postsynId], ...
    [synapses.syn_edgeIdx], edges);

synapses = struct2table(synapses);
if nargout > 1
    boutons = cat(2, boutons, synapses);
end

end

function presynId = getPresynId(postsynId, edgeIdx, edges)

presynId = cellfun(@(idx, ps) setdiff(edges(idx,:), ps), ...
    vertcat(edgeIdx{:}), vertcat(postsynId{:}), 'uni', 0);
presynId = mat2cell(presynId, cellfun(@length, edgeIdx), 1);

end

