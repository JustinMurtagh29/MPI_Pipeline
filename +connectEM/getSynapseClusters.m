function [nBoutons,meanSynPerBouton] = getSynapseClusters(p,agglos,show,clusterDistance)
% p = Gaba.getSegParameters('ex145_ROI2017');
% [graph, segmentMeta, borderMeta] = Seg.IO.loadGraph(p, false);
if ~exist('show','var') || isempty(show)
    show = 0;
end
if ~exist('clusterDistance','var') || isempty(clusterDistance)
    clusterDistance = 1000; %nm
end
load(fullfile(p.saveFolder,'connectomeState','SynapseAgglos_v2.mat'),'synapses') % (erzeugt via E:\workspace\pipeline\Benedikt\+L4\+Synapses\+Scripts\synapseDetection.m)

 
presynSegIds = cellfun(@(x) x(1),synapses.presynId); % only one presyn id necessary
nBoutons = NaN(numel(agglos),1);
meanSynPerBouton = nBoutons;
for s = 1:numel(agglos)
    [~,ind] = ismember(presynSegIds,agglos(s).nodes(:,4));
    presynSegCoords = bsxfun(@times,agglos(s).nodes(ind(ind~=0),1:3),[11.24,11.24,28]);
    if size(presynSegCoords,1) > 1
        synClusters = clusterdata(presynSegCoords, ...
            'linkage', 'single', 'criterion', 'distance', 'cutoff', clusterDistance);
    else
        synClusters = 1;
    end
    nBoutons(s) = max(synClusters);
    meanSynPerBouton(s) = mean(histc(synClusters,1:nBoutons(s)));
    if show
        figure; hold all
        skel = Superagglos.toSkel(agglos(s));
        skel.plot;
        scatter3(presynSegCoords(:,1),presynSegCoords(:,2),presynSegCoords(:,3),150,synClusters,'o','filled')
        colormap lines
    end
end


% Dafür wurden schon einzelne Interfaces zu Synapsen geclustered, wenn diese nach einer lokalen Agglomeration den gleichen prä- und postsynaptischen Partner haben.
% 
% Die edgeIdx der Synapsen beziehen sich auf die edges in graph.mat.
% Falls du die einzelnen Interfaces laden willst, dann von der Datei globalSynScores.mat. Threshold -1.67. Diese beziehen sich auf die edges in globalEdges.mat.

 