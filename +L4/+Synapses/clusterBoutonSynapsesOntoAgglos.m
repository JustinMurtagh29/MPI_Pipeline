function synOntoPS = clusterBoutonSynapsesOntoAgglos( boutons, agglosPS, ...
    borderCom, d, b2PS, d_nPS)
%CLUSTERBOUTONSYNAPSESONTOAGGLOS Cluster the synapses of a bouton and a
% postsynaptic agglo.
% INPUT boutons: table
%           Table of boutons containing the rows edgeIdx and postsynId.
%       agglosPS: [Nx1] cell
%           Cell array of postsynaptic agglos.
%       borderCom: [Nx3] single
%           Global border center of mass.
%       d: double
%           Distance below which points get clustered.
%       b2PS: [Nx1] cell
%           Mapping of boutons onto postsynaptic target.
%       d_nPS: (Optional) double
%           Agglomeration distance for synapses that are not onto
%           postsynaptic agglos.
%           (Default: d)
% OUTPUT synOntoPS: struct
%           Struct containing the fields 'n_syn_onto_PS', 'edgeIdx', and
%           'targetIdx' that contains the number of synapses, the global
%           edge index for all interfaces of all synapse and the index of
%           all target agglos of all synapses for each bouton.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('d_nPS', 'var') || isempty(d_nPS)
    d_nPS = d;
end

if isinteger(borderCom)
    borderCom = single(borderCom);
end

m = max([max(cellfun(@max, boutons.postsynId)), ...
         max(cellfun(@max, agglosPS))]);
idToAgglo = Seg.Global.eClassLookup(agglosPS, m, false);

idx = find(~cellfun(@isempty, b2PS));
synOntoPS.n_syn_onto_PS = zeros(length(b2PS), 1);
synOntoPS.n_syn = zeros(length(b2PS), 1);
synOntoPS.edgeIdx = cell(length(b2PS), 1);
synOntoPS.targetIdx = cell(length(b2PS), 1);
for i = 1:length(idx)
    psIdx = b2PS{idx(i)};
    group = zeros(length(boutons.postsynId{idx(i)}), 1);
    target = zeros(size(group));
    % for each postsynaptic target
    for j = 1:length(psIdx)
        isToPS = ismember(idToAgglo(boutons.postsynId{idx(i)}), psIdx(j));
        target(isToPS) = psIdx(j);
        if sum(isToPS) > 1
            T = clusterdata(borderCom(boutons.edgeIdx{idx(i)}(isToPS), :), ...
                'linkage', 'single', 'criterion', 'distance', ...
                'cutoff', d);
            group(isToPS) = max(group) + T;
        else
            group(isToPS) = max(group) + 1;
        end
    end
    
    n_synOntoPS = max(group);
    
    % cluster remaining synapses
    isToPS = ~ismember(idToAgglo(boutons.postsynId{idx(i)}), psIdx);
    target(isToPS) = 0;
    if sum(isToPS) > 1
        T = clusterdata(borderCom(boutons.edgeIdx{idx(i)}(isToPS), :), ...
            'linkage', 'single', 'criterion', 'distance', ...
            'cutoff', d_nPS);
        group(isToPS) = max(group) + T;
    else
        group(isToPS) = max(group) + 1;
    end
    
    synInt = boutons.edgeIdx{idx(i)}(group > 0);
    synOntoPS.n_syn(idx(i)) = max(group);
    synOntoPS.n_syn_onto_PS(idx(i)) = n_synOntoPS;
    synOntoPS.edgeIdx{idx(i)} = accumarray(group(group > 0), ...
        synInt, [], @(x){x});
    synOntoPS.targetIdx{idx(i)} = accumarray(group(group > 0), ...
        target(group > 0), [], @(x)x(1));
end


end

