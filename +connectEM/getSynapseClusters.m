function [nClusters,meanSynPerCluster,meanClusterVolumes] = getSynapseClusters(p,agglos,synState,show,clusterDistance)

if ~exist('synState','var') || isempty(synState)
    synState = 'SynapseAgglos_v2';
end

if ~exist('show','var') || isempty(show)
    show = 0;
end
if ~exist('clusterDistance','var') || isempty(clusterDistance)
    clusterDistance = 1000; %nm
end
load(fullfile(p.saveFolder,'connectomeState',[synState,'.mat']),'synapses') % (erzeugt via E:\workspace\pipeline\Benedikt\+L4\+Synapses\+Scripts\synapseDetection.m)
segmentMeta = load(fullfile(p.saveFolder, 'segmentMeta.mat'),'voxelCount');

 
presynSegIds = cellfun(@(x) x(1),synapses.presynId); % only one presyn id necessary
presynSegVol = cellfun(@(x) sum(segmentMeta.voxelCount(x)),synapses.presynId); % get vol of presyn
nClusters = NaN(numel(agglos),1);
meanSynPerCluster = nClusters;
meanClusterVolumes = nClusters;
for s = 1:numel(agglos)
    [presynInd,ind] = ismember(presynSegIds,agglos(s).nodes(:,4));
    presynSegCoords = bsxfun(@times,agglos(s).nodes(ind(ind~=0),1:3),[11.24,11.24,28]);
    if size(presynSegCoords,1) > 1
        synClusters = clusterdata(presynSegCoords, ...
            'linkage', 'single', 'criterion', 'distance', 'cutoff', clusterDistance);
    elseif size(presynSegCoords,1)==1
        synClusters = 1;
    else
        nClusters(s) = 0;
        meanClusterVolumes(s) = 0;
        meanSynPerCluster(s) = 0;
        continue
    end
    meanClusterVolumes(s) = mean(accumarray(synClusters,presynSegVol(presynInd))*prod([0.01124,0.01124,0.028]));
    nClusters(s) = max(synClusters);
    meanSynPerCluster(s) = mean(histc(synClusters,1:nClusters(s)));
    if show
        figure; hold all
        skel = Superagglos.toSkel(agglos(s));
        skel.plot;
        scatter3(presynSegCoords(:,1),presynSegCoords(:,2),presynSegCoords(:,3),150,synClusters,'o','filled')
        colormap lines
    end
end