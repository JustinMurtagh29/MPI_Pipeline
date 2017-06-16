load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
load('/gaba/scratch/mberning/20170404T172554_agglomerationNew/aggloFinal.mat', 'final');
m = load([p.saveFolder 'globalEdges.mat']);
edges = m.edges;
m = load([p.saveFolder 'globalSynScores.mat']);
synScores = nan(size(edges, 1), 2);
synScores(m.edgeIdx,:) = m.synScores;
denEdgeIdx = any(ismember(edges, final.dendrites{7}), 2);
denEdges = edges(denEdgeIdx,:);
denScores = synScores(denEdgeIdx, :);
denSyn = any(denScores > -1.23, 2);
m = load([p.saveFolder 'globalBorder.mat']);
borderCoM = m.borderCoM;
clear m
edgeIdx = ismember(edges, denEdges(denSyn,:), 'rows');
theseSynapticLocations = borderCoM(edgeIdx,:);
% Single linkage clustering
dist = squareform(pdist(bsxfun(@times, double(theseSynapticLocations), p.raw.voxelSize)));
[row, col] = find(dist < 350);
cc = Graph.findConnectedComponents([row col], false, false);
theseSynapticLocationsClustered = cellfun(@(x)mean(theseSynapticLocations(x,:),1), cc, 'uni', 0);
theseSynapticLocationsClustered = cell2mat(theseSynapticLocationsClustered);
fid = fopen('/home/mberning/Desktop/synapseLocations.txt', 'w');
for i=1:size(theseSynapticLocationsClustered,1)
    fprintf(fid, '%i, %i, %i\n', round(theseSynapticLocationsClustered(i,:)));
end
fclose(fid);
%save('/home/mberning/localStorage/forConnectome.mat');

%% Connectome
synapticEdges = edges(any(synScores > -1.23,2),:);
synapticCoM = borderCoM(any(synScores > -1.23,2),:);
% Create matrix of size of synaptic edges to 
linkingComponentsD = zeros(size(synapticEdges));
nrDendrites = length(final.dendrites);
for i=1:nrDendrites
    linkingComponentsD(ismember(synapticEdges, final.dendrites{i})) = i;
end
linkingComponentsA = zeros(size(synapticEdges));
nrAxons = length(final.axons);
for i=1:nrAxons
    linkingComponentsA(ismember(synapticEdges, final.axons{i})) = i;
end
% Remove synapses detected within dendrites (1%) and axons (0.1%)
linkingComponentsD(all(linkingComponentsD > 0, 2),:) = 0;
linkingComponentsA(all(linkingComponentsA > 0, 2),:) = 0;
% Keep only synapses that connect current state of axon and dendrites
componentLinks = cat(2, linkingComponentsD, linkingComponentsA);
idx = sum(componentLinks > 0, 2) == 2;
theseComponentLinks = componentLinks(idx,:);
% Single linkage clustering
dist = squareform(pdist(bsxfun(@times, double(synapticCoM(idx,:)), p.raw.voxelSize)));
[row, col] = find(dist < 350);
cc = Graph.findConnectedComponents([row col], false, false);
% Test whether CC usually link the same component (checked a few and looks
% OK, found no contradictory links, investigation should probably be
% increased
ccSize = cellfun(@numel, cc);
find(ccSize > 1, 'first)')
theseComponentLinks(cc{12},:)
% Keep only one synapse of each CC
ccRep = cellfun(@(x)x(1), cc);
theseComponentLinks = theseComponentLinks(ccRep,:);
% Sort out zeros, should be a better way to do this
theseComponentLinksNew = zeros(size(theseComponentLinks,1),2);
for i=1:size(theseComponentLinks,1)
    theseComponentLinksNew(i,:) = theseComponentLinks(i,theseComponentLinks(i,:) > 0);
end
connectome = accumarray(theseComponentLinksNew, 1, ...
    [numel(final.dendrites) numel(final.axons)], @sum);
%save('/home/mberning/localStorage/forConnectome2.mat');

%% Try to plot
figure;
[x, y, v] = find(connectome(randperm(size(connectome, 1)),randperm(size(connectome, 2))));
v(v > 4) = 4;
colors = [0.7 0.7 0.7; 0.5 0.5 0.5; 0.5 0.2 0.6; 0 0 0];
scatterhist(x,y,'Group',v,'Location','SouthEast',...
    'Direction','out','Color',colors,'LineStyle',{'-', ':', '-.', '--'},...
    'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[3,3,3,3], ...
    'NBins', size(connectome));
xlabel('dendrites');
ylabel('axons');
xlim([1 size(connectome, 1)]);
ylim([1 size(connectome, 2)]);
legend('1', '2', '3', '>4');
axis off;

[x, y, v] = find(connectome);
v(v > 4) = 4;
colors = [0.7 0.7 0.7; 0.5 0.5 0.5; 0.5 0.2 0.6; 0 0 0];
figure;
scatterhist(x,y,'Group',v,'Location','SouthEast',...
    'Direction','out','Color',colors,'LineStyle',{'-', ':', '-.', '--'},...
    'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[3,3,3,3], ...
    'NBins', size(connectome));
xlabel('dendrites');
ylabel('axons');
xlim([1 size(connectome, 1)]);
ylim([1 size(connectome, 2)]);
legend('1', '2', '3', '>4');
axis off;

%% Plot with black background

figure;
[x, y, v] = find(connectome(randperm(size(connectome, 1)),randperm(size(connectome, 2))));
v(v > 4) = 4;
colors = [0.9 0.5 0.5; 0.5 0.9 0.5; 0.5 0.5 0.9; 1 1 1];
figure;
scatterhist(x,y,'Group',v,'Location','SouthEast',...
    'Direction','out','Color',colors,'LineStyle',{'-', ':', '-.', '--'},...
    'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[3,3,3,3], ...
    'NBins', round(size(connectome)./100));
axis off;
set(gcf, 'Color', [0 0 0]);
