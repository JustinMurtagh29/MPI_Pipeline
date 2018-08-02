% proofreading of whole cell axons (2018-08-02)
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo();


%% load agglos

p = Gaba.getSegParameters('ex145_ROI2017');
wcFile = fullfile(p.agglo.saveFolder, ...
    'wholeCells_GTAxon_08_v4_splitHealed_v1.mat');
somaFile = fullfile(p.agglo.saveFolder, 'somata_07.mat');

m = load(wcFile);
wcAgglos = m.wholeCells;
wcAgglosC = SuperAgglo.clean(wcAgglos);

m = load(somaFile);
somaAgglos = m.somata;


%% get somata and axons

soma = connectEM.WholeCell.getSoma(wcAgglosC, somaAgglos);


%% write to WK

% somata
skel = Superagglos.toSkel(Superagglos.deleteNodes(wcAgglosC, soma, true));
skel.names = arrayfun(@(x)sprintf('WC_%02d_Soma', 1:length(wcAgglosC), ...
    'uni', 0));

% axons
skel = Superagglos.toSkel(Superagglos.deleteNodes(wcAgglosC, ...
    {wcAgglosC.axons}, true), skel);
skel.names = arrayfun(@(x)sprintf('WC_%02d_Axon', 1:length(wcAgglosC), ...
    'uni', 0));

