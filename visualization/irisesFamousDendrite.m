spineApp = parseNml('C:\Users\mberning\Downloads\DendriteComm.nml');
synTyp = parseNml('C:\Users\mberning\Downloads\SynTypComm.nml');

%% Find dendrite in syn typ file
for i = 1:length(synTyp)
    display([num2str(i) ': ' synTyp{i}.name]);
end
dendrite = synTyp(18);

%% Find order of nodes
uniqueIdsInEdges = unique(dendrite{1}.edges(:)); % only necessary if nodes IDs are missing
idCountInEdges = histc( dendrite{1}.edges(:), uniqueIdsInEdges);
idxOrder3 = find(idCountInEdges == 3);
idxOrder1 = find(idCountInEdges == 1);

%% Look at it
hist(idCountInEdges); % there are 14 of order 4 and 1 of order 5, what to do with that? currenlty ignored

%% or that
figure; hold on;
for i=1:length(dendrite{1}.edges)
    X = [dendrite{1}.nodes(dendrite{1}.edges(i,1),1) dendrite{1}.nodes(dendrite{1}.edges(i,2),1)];
    Y = [dendrite{1}.nodes(dendrite{1}.edges(i,1),2) dendrite{1}.nodes(dendrite{1}.edges(i,2),2)];
    Z = [dendrite{1}.nodes(dendrite{1}.edges(i,1),3) dendrite{1}.nodes(dendrite{1}.edges(i,2),3)];
    plot3(X,Y,Z, 'LineWidth', 1);
end

%% find position of all nodes of order 1 (to use as reference for path length, taken from plot above because it is one end probably in soma)
startIdx = find(all(dendrite{1}.nodes(:,1:3) == repmat([2649 2409 2121],size(dendrite{1}.nodes,1),1),2));

%% constuct sparse distance matrix
addpath('C:\code\gaimc\'); % For dijkstra: http://www.mathworks.com/matlabcentral/fileexchange/24134-gaimc-graph-algorithms-in-matlab-code
A = sparse(size(dendrite{1}.nodes,1));
for j=1:size(dendrite{1}.edges,1)
    id1 = dendrite{1}.edges(j,1);
    id2 = dendrite{1}.edges(j,2);
    dist = norm((dendrite{1}.nodes(id1,1:3) - dendrite{1}.nodes(id2,1:3)).*[11.28 11.28 28]);
    A(id1, id2) = dist;
    A(id2, id1) = dist;
end
% Find distances to soma node
[distances preVec] = dijkstra(A, startIdx);

%% Iterate over all nodes to find comments and collect type and linearized distance to soma
idx = 1;
for i=1:length(dendrite{1}.nodesAsStruct)
    if ~isempty(dendrite{1}.nodesAsStruct{i}.comment)
        foundSth(idx).comment = dendrite{1}.nodesAsStruct{i}.comment;
        foundSth(idx).nodeID = dendrite{1}.nodesAsStruct{i}.id;
        foundSth(idx).nodeOrder = idCountInEdges(uniqueIdsInEdges(i));
        foundSth(idx).distanceToSoma = distances(uniqueIdsInEdges(i));
        % weird id in nodeStruct is already class char
        display(['Found comment: ' foundSth(idx).comment ' at node id ' foundSth(idx).nodeID ' ( node order ' num2str(foundSth(idx).nodeOrder) ')']);      
        % aftertought: find indicator for synapse type
        foundSth(idx).synTyp = dendrite{1}.nodesAsStruct{i}.comment(end-1:end);
        foundSth(idx).preSynTree = str2num(dendrite{1}.nodesAsStruct{i}.comment(2:3));
        idx = idx + 1;
    end
end

%% sort out comments that do not indicate synapse type

cc = []; tc = []; inh = [];
for i=1:length(foundSth)
    if all(~strcmp(foundSth(i).synTyp,{'CC' 'TC' 'NH'}));
        display(['Ignoring comment: ' foundSth(i).synTyp ' at node id ' foundSth(i).nodeID]);
    else
        switch foundSth(i).synTyp
            case 'CC'
                cc(end+1,1) = foundSth(i).distanceToSoma;
                cc(end,2) = foundSth(i).preSynTree;
            case 'TC'
                tc(end+1,1) = foundSth(i).distanceToSoma;
                tc(end,2) = foundSth(i).preSynTree;
            case 'NH'
                inh(end+1,1) = foundSth(i).distanceToSoma;
                inh(end,2) = foundSth(i).preSynTree;
        end
    end
end

%% a first try at visualization (maybe needs redoing after skeleton (or comment out second row of cc,tc,inh)
figure; hold on;
nrBins = 20;
bins = min(distances):(max(distances)-min(distances))/nrBins:max(distances);
a = hist(cc,bins);
b = hist(tc,bins);
c = hist(inh,bins);
bar(bins, [a; b; c]', 'stacked');
xlim([min(distances) max(distances)]);
legend('Cortiocortical', 'Thalamocortical', 'Inhibitory');
xlabel('distance along dendrite measured from node (i guess center) of soma');
ylabel('occurences');

%% Save to somewhere
save('C:\Users\mberning\Desktop\SynapseDistirbutionAlondDendrite.mat');

%% Save all dendrites to 3 files according to type
clc;
for i = 1:length(synTyp)
    if strcmp(synTyp{i}.name(1:4), 'Tree') 
        treeComment(i) = str2num(synTyp{i}.name(end-1:end));
        display([num2str(i) ': ' synTyp{i}.name]);
    end
end

ccSkel = {};
for i=1:size(cc,1)
    idx = find(treeComment == cc(i,2));
    ccSkel(end+1) = synTyp(idx);
end

tcSkel = {};
for i=1:size(tc,1)
    idx = find(treeComment == tc(i,2));
    tcSkel(end+1) = synTyp(idx);
end

inhSkel = {};
for i=1:size(inh,1)
    idx = find(treeComment == inh(i,2));
    inhSkel(end+1) = synTyp(idx);
end

convertKnossosNmlToHoc2(ccSkel, ['I:\CortexConnectomics\Manuel\sync\wholeCell\DendriteFromIris\cc.hoc'], 0, 1, 0, 0, [11.28 11.28 28]);
convertKnossosNmlToHoc2(tcSkel, ['I:\CortexConnectomics\Manuel\sync\wholeCell\DendriteFromIris\tc.hoc'], 0, 1, 0, 0, [11.28 11.28 28]);
convertKnossosNmlToHoc2(inhSkel, ['I:\CortexConnectomics\Manuel\sync\wholeCell\DendriteFromIris\inh.hoc'], 0, 1, 0, 0, [11.28 11.28 28]);
