% script to determine apical dendrite specificity
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

plotting = false; %#ok<*UNRCH>

%% load apicals, axons and boutons

% apicals
m = load('/tmpscratch/zecevicm/L4/apicalDendrites/apicalDendritesL4.mat');
apicalDenAgglos = m.filteredAgglos;

% axons
m = load('/tmpscratch/mberning/axonQueryResults/postQueryAgglomerates.mat');
axonAgglos.segIds = m.axonsPostQuery(m.isAbove5um);
axonAgglos.segIds  = axonAgglos.segIds(2:end); % ignore super-merger

% boutons
m = load('/gaba/u/bstaffle/data/L4/Synapses/boutonAgglos_v3.mat');
boutons = m.boutons_nOv;

%% get boutons onto apicals

[bIdx, b2Apical] = L4.Synapses.getBoutonsOntoAgglo(boutons.postsynId, ...
    apicalDenAgglos);

% list which bouton is onto which apical
isOntoApical = zeros(size(boutons, 1), 1);
isOntoApical(cell2mat(bIdx)) = repelem(1:length(bIdx), cellfun(@length, bIdx));

%% cluster synapses

load('/gaba/u/bstaffle/code/workspace/allParameter20170217.mat')
m = load(fullfile(p.saveFolder, 'globalBorder.mat'), 'borderCoM');
borderCom = m.borderCoM;

synapsesOntoApicals = L4.Synapses.clusterBoutonSynapsesOntoAgglos( ...
    boutons, apicalDenAgglos, borderCom, 30, b2Apical);
%% get boutons for axons

[b2agglo, axonAgglos.boutons] = L4.Synapses.getCollectedBoutons( ...
    axonAgglos.segIds , boutons.agglo );

%% extend boutons table
boutons = cat(2, boutons, table(b2agglo, 'VariableNames', {'AxonAggloIdx'}));
boutons = cat(2, boutons, table(synapsesOntoApicals.targetIdx, ...
    synapsesOntoApicals.n_syn_onto_PS, ...
    'VariableNames', {'PostsynIdx', 'n_syn_onto_PS'}));
axonAgglos.boutonsOnApical = cellfun(@(x)isOntoApical(x), ...
    axonAgglos.boutons, 'uni', 0);
fractionSynapsesOnApical = cellfun( ...
    @(x)(sum(synapsesOntoApicals.n_syn_onto_PS(x)) - 1)/ ...
         (sum(synapsesOntoApicals.n_syn(x)) - 1), ...
         axonAgglos.boutons(axonAgglos.hasBoutonsOntoApicals));
synapsesOnApical = cellfun( ...
    @(x)(sum(synapsesOntoApicals.n_syn_onto_PS(x)) - 1), ...
         axonAgglos.boutons(axonAgglos.hasBoutonsOntoApicals));
axonAgglos.hasBoutonsOntoApicals = cellfun(@any, axonAgglos.boutonsOnApical);
axonAgglos.fractionOntoApicals = cellfun(@(x)sum(x>0)/length(x), ...
    axonAgglos.boutonsOnApical);
fractionOntoApicalMinusOne = cellfun(@(x)(sum(x>0) - 1)/(length(x) - 1), ...
    axonAgglos.boutonsOnApical(axonAgglos.hasBoutonsOntoApicals));
totalNumber = cellfun(@(x)length(x) - 1, ...
    axonAgglos.boutonsOnApical(axonAgglos.hasBoutonsOntoApicals));

%% fractional specificity

if plotting 
    figure
    histogram(fractionOntoApicalMinusOne)
    xlabel('Fraction of synapses onto apicals (minus one)')
    ylabel('Count')
    Visualization.Figure.plotDefaultSettings()
    a = gca;
    a.XLim = [0, 1];

    figure
    scatter(fractionOntoApicalMinusOne, totalNumber);
    xlabel('Fraction of synapses onto apicals (minus one)')
    ylabel('Total number of synapses (minus one)')
    Visualization.Figure.plotDefaultSettings()
    a = gca;
    a.XLim = [0, 1];
    a.YScale = 'log';
end


%% example to nml

if plotting
    load('/gaba/u/bstaffle/code/workspace/allParameter20170217.mat')
    load([p.saveFolder 'segmentMeta.mat'], 'point')
    point = point';
    skel = L4.Agglo.agglo2Nml(cat(1, apicalDenAgglos(1), ...
        boutons.agglo(bIdx{1})), point);
    skel.names{1} = 'Apical';
    for i = 2:skel.numTrees()
        skel.names{i} = sprintf('Bouton_%03d', i - 1);
    end
    skel.write('BoutonsOntoApicalExample.nml')
end
