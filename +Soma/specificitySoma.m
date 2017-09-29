% script to determine soma dendrite specificity (using synapses)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>, modified by Robin Hesse

plotting = false; %#ok<*UNRCH>


%% load somas, axons and synapses

load(fullfile('/home/robin/gaba/gaba/u/mberning/results/pipeline/20170217_ROI', 'allParameter.mat'));

% somas
m = load('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/somas_with_merged_somas.mat');
somaAgglos = m.somas(:,3);

% axons
m = load('/gaba/tmpscratch/mberning/axonQueryResults/postQueryAgglomerates.mat');
axonAgglos.segIds = m.axonsPostQuery(m.isAbove5um);
axonAgglos.segIds  = axonAgglos.segIds(2:end); % ignore super-merger

% synapses
m = load('/gaba/u/bstaffle/data/L4/Synapses/SynapseAgglos_v1.mat');
synapses = m.synapses;


%% get synapse agglo mappings

[syn2Agglo, pre2syn, post2syn, connectome] = ...
    L4.Synapses.synapsesAggloMapping(synapses, axonAgglos.segIds, ...
    somaAgglos);


%% examples to nml (connections with > 1 synapses)

if plotting
    p = Gaba.getSegParameters('ex145_ROI2017');
    [graph, segmentMeta, borderMeta] = Seg.IO.loadGraph(p, false);
    
    idx = find(cellfun(@length, connectome.synIdx) > 1);
    i = 4;
    skel = L4.Agglo.agglo2Nml( ...
        [axonAgglos.segIds(connectome.edges(idx(i), 1)); ...
        somaAgglos(connectome.edges(idx(i), 2))], segmentMeta.point);
    skel = L4.Agglo.agglo2Nml( ...
        synapses.edgeIdx(connectome.synIdx{idx(i)}, 1), ...
        borderMeta.borderCoM, skel);
    skel.write('ConnectoneEntrySample.nml')
end

%% specificity

% axons onto somas
axonIdx = unique(connectome.edges(:,1));

% get synapses onto somas
isSynOntoSoma = false(size(synapses, 1), 1);
isSynOntoSoma(cell2mat(connectome.synIdx)) = true;

% get synapse counts for axons
totalNoSyn = cellfun(@length, pre2syn(axonIdx));
synOntoSomas = cellfun(@(x)sum(isSynOntoSoma(x)), pre2syn(axonIdx));

% specificities
sp_minus_one = (synOntoSomas - 1)./(totalNoSyn - 1);

% all information to table
axonsOntoSomas = table(axonIdx, totalNoSyn, synOntoSomas, ...
    sp_minus_one, 'VariableNames', {'axonIdx', 'totalNoSyn', ...
    'synOntoSomas', 'specificity_minus_one'});

% remopve all axons that have only one synapse & calculate avg specificity
axonsOntoSomasWithoutAxonsThatHaveOnlyOneSynapse = ...
    axonsOntoSomas(axonsOntoSomas.totalNoSyn~=1,:);

avgSp = mean(axonsOntoSomasWithoutAxonsThatHaveOnlyOneSynapse.specificity_minus_one);
s = std(axonsOntoSomasWithoutAxonsThatHaveOnlyOneSynapse.specificity_minus_one);

%% fractional specificity

if plotting 
    figure
    h = histogram(axonsOntoSomas.specificity_minus_one)
    xlabel('Fraction of synapses onto somas (minus one)')
    ylabel('Count')
    Visualization.Figure.plotDefaultSettings()
    a = gca;
    a.XLim = [0, 1];
    
    
    hold on 
    
    x = 0:0.05:1;
    %y(1,1:21) = avgSp;
    y(1,1:21) = mean(h.Values);
    %err(1,1:21) = s;
    err(1,1:21) = std(h.Values);
    errorbar(x,y,err)

    hold off
    
    
    figure
    h = histogram(axonsOntoSomas.specificity_minus_one)
    xlabel('Fraction of synapses onto somas (minus one)')
    ylabel('Count')
    Visualization.Figure.plotDefaultSettings()
    a = gca;
    a.XLim = [0, 1];
    a.YScale = 'log';
    
    
    hold on 
    
    x = 0:0.05:1;
    y(1,1:21) = avgSp;
    %y(1,1:21) = mean(h.Values);
    %err(1,1:21) = s;   %makes no sense in log scale
    %err(1,1:21) = std(h.Values);
    plot(x,y)

    hold off
    
    x = 0:0.05:1;
    y(1,1:21) = avgSp;
    %y(1,1:21) = mean(h.Values);
    err(1,1:21) = s;   %makes no sense in log scale
    %err(1,1:21) = std(h.Values);
    errorbar(x,y,err)
    
    figure
    scatter(axonsOntoSomas.specificity_minus_one, ...
        axonsOntoSomas.totalNoSyn);
    xlabel('Fraction of synapses onto somas (minus one)')
    ylabel('Total number of synapses (minus one)')
    Visualization.Figure.plotDefaultSettings()
    a = gca;
    a.XLim = [0, 1];
    a.YScale = 'log';
end

