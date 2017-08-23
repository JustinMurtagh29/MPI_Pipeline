function synapses = clusterBoutonSynapses( boutons, borderCom, d, ...
    voxelSize )
%CLUSTERBOUTONSYNAPSES Cluster synapses of a bouton based on distance.
% INPUT boutons: table
%           Table of boutons containing the rows edgeIdx and postsynId.
%       borderCom: [Nx3] single
%           Global border center of mass.
%       d: double
%           Distance below which points get clustered.
%       voxelSize: (Optional) [1x3] double
%           The voxel size.
%           (Default: [1, 1, 1])
% OUTPUT synapses: struct
%           Struct containing the fields 'n_syn', 'edgeIdx', and
%           'targetIdx' that contains the number of synapses, the global
%           edge index for all interfaces of all synapse and the index of
%           all target agglos of all synapses for each bouton.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if isinteger(borderCom)
    borderCom = single(borderCom);
end

if exist('voxelSize', 'var') && ~isempty(voxelSize)
    borderCom = bsxfun(@times, borderCom, voxelSize(:)');
end

noB = size(boutons, 1);
synapses.n_syn = zeros(noB, 1);
synapses.edgeIdx = cell(noB, 1);
synapses.targetIdx = cell(noB, 1);
for i = 1:noB
    edgeIdx = boutons.edgeIdx{i};
    if length(edgeIdx) > 1
        T = clusterdata(borderCom(edgeIdx, :), ...
            'linkage', 'single', 'criterion', 'distance', 'cutoff', d);
    else
        T = 1;
    end
    synapses.n_syn(i) = max(T);
    synapses.edgeIdx{i} = accumarray(T, edgeIdx, [], @(x){x});
    synapses.targetIdx{i} = accumarray(T, ...
        boutons.postsynId{i}, [], @(x)x(1));
end


end

