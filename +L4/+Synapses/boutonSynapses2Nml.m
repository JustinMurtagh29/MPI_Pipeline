function skel = boutonSynapses2Nml( synapses, borderCom )
%BOUTONSYNAPSES2NML Write synapses of boutons to nml.
% INPUT synapses: [Nx1] cell of [Mx1] cell
%           Cell array where each entry contains all synapses for one
%           bouton. Each synapse is another cell array containing the
%           synaptic interfaces indices of the synapse.
%       borderCom: [Nx3] double
%           The global edge com list.
% OUTPUT skel: skeleton object
%           Skeleton object containing each agglomeration as a single tree.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

skel = Skeleton.fromMST(cellfun(@(ids) {borderCom(ids, :)}, ...
    vertcat(synapses{:})), [11.24, 11.24, 28]);
skel = Skeleton.setParams4Dataset(skel, 'ex145_ROI2017');
bIdx = repelem(1:length(synapses), cellfun(@length, synapses))';
skel.names = arrayfun(@(x)sprintf('Bouton%03d', x), bIdx, 'uni', 0);
cmap = lines();
cmap(:,4) = 1;
skel.colors = num2cell(cmap(bIdx,:), 2);


end

