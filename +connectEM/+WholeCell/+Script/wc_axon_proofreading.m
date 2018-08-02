% proofreading of whole cell axons (2018-08-02)
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo();


%% load agglos

p = Gaba.getSegParameters('ex145_ROI2017');
wcFile = fullfile(p.agglo.saveFolder, ...
    'wholeCells_GTAxon_08_v5.mat');
somaFile = fullfile(p.agglo.saveFolder, 'somata_07.mat');

m = load(wcFile);
wcAgglos = m.wholeCells;
wcAgglosC = SuperAgglo.clean(wcAgglos);
numWC = length(wcAgglos);

m = load(somaFile);
somaAgglos = m.somata;


%% get somata and axons

soma = connectEM.WholeCell.getSoma(wcAgglosC, somaAgglos);
for i = 1:length(wcAgglosC)
    if ~islogical(wcAgglosC(i).axon)
        if all(isnan(wcAgglosC(i).axon))
            wcAgglosC(i).axon = false(size(wcAgglosC(i).nodes, 1), 1);
        else
            error('Axons that are not fully nan.');
        end
    end
end


%% write to WK

% somata
wcSomas = SuperAgglo.toMST( ...
    Superagglos.deleteNodes(wcAgglosC, soma, true), [11.24, 11.24, 28]);
skel = Superagglos.toSkel(wcSomas);
skel.names(1:numWC) = arrayfun(@(x)sprintf('WC_%02d_Soma', x), ...
    1:length(wcAgglosC), 'uni', 0);

% axons
wcAxons = SuperAgglo.toMST( ...
    Superagglos.deleteNodes(wcAgglosC, {wcAgglosC.axon}, true), ...
    [11.24, 11.24, 28]);
skel = Superagglos.toSkel(wcAxons, skel);
skel.names(numWC + 1:end) = arrayfun(@(x)sprintf('WC_%02d_Axon', x), ...
    1:length(wcAgglosC), 'uni', 0);

[~, wcName] = fileparts(wcFile);

skel = skel.setDescription( ...
    sprintf(['Whole cell (%s) axons and somata. ' ...
    '(function: %s, git hash: %s)'], wcName, info.filename, ...
    info.git_repos{1}.hash));

% sort by wc
idx = reshape(1:skel.numTrees(), 2, [])';
skel.thingIDs = idx(:);
skel = skel.sortTreesById();
skel.colors(1:2:end) = {[0 0 1 1]};
skel.colors(2:2:end) = {[1 0 0 1]};

for i = 1:numWC
    [skel, gid] = skel.addGroup(sprintf('WC_%02d', i));
    skel = skel.addTreesToGroup((2*i-1):(2*i), gid);
end

