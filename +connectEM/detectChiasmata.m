function output = detectChiasmata(p, nodesV, edges, visualize, outputFolder )
% Detect chiasmata in skeletons based on marching sphere algorithm
% Nodes should be in voxel, scaled here
%load('/gaba/u/kboerg/biggest/biggestpre1.mat');
% Create output folder if it does not exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Scale to nm
% Make sure edges are unique
[edges, ~, idxE] = unique(edges, 'rows');

% Some parameter for algorithm
p.sphereRadiusOuter = 10000; % in nm
p.sphereRadiusInner = 1000; % in nm
p.voxelSize = [11.24 11.24 28];
nodes=bsxfun(@times,nodesV,p.voxelSize);
% for each node ("marching sphere" approach to merger detection)
isIntersection = false(size(nodes,1),1);
nrExits = zeros(size(nodes,1),1);
if size(nodes, 1) < 1E6
    for i=1:size(nodes,1)
        [isIntersection(i),nrExits(i)] = detectChiasmataSub(i,nodes,edges,p);
    end
else
    parfor i=1:size(nodes,1)
        [isIntersection(i),nrExits(i)] = detectChiasmataSub(i,nodes,edges,p);
    end
end
%save(['/gaba/u/kboerg/biggest/2biggest' num2str(startidx) '.mat'], 'isIntersection');

% Find CC of detected intersections according to graph
[output, queryIdx] = connectEM.detectChiasmataPostSingleNodeLabel(edges, isIntersection, nrExits, nodes, p, nodesV, ones(size(edges,1),1));

if visualize
    % Write result to skletons for control (detection of intersections)
    comments = cell(size(nodesV,1), 1);
    comments(isIntersection) = arrayfun(@(x)strcat('intersection, ', num2str(x), ' exits'), nrExits(isIntersection), 'uni', 0);
    writeNml([outputFolder 'skelDetected.nml'], writeSkeletonFromNodesAndEdges({nodesV}, {edges}, {comments}, {'axon'}, {[0 0 1 1]}));
    % Write result to skletons for control (query visualization)
    clear comments;
    n{1} = nodesV(~isIntersection,:);
    % Keep only nodes & edges not part of any intersection
    edgeOffset = cumsum(accumarray(find(~isIntersection), 1))';
    e{1} = edges(~any(ismember(edges, find(isIntersection)),2),:);
    e{1} = edgeOffset(e{1});
    comments{1} = cell(size(n{1},1), 1);
    t{1} = 'original tree without detected intersections';
    col{1} = [0 0 1 1];
    for idx=1:length(queryIdx)
        n{idx+1} = cat(1, nodesV(output.ccCenterIdx(idx),:), nodesV(output.queryIdx{idx},:));
        e{idx+1} = cat(2, ones(size(n{idx+1},1)-1,1), (2:size(n{idx+1},1))');
        comments{idx+1} = cell(size(n{idx+1},1), 1);
        comments{idx+1}{1} = 'center node of intersection, edges inidcating queries';
        t{idx+1} = ['intersection ' num2str(idx)];
        col{idx+1} = [1 0 0 1];
    end
    writeNml([outputFolder 'skelQueries.nml'], writeSkeletonFromNodesAndEdges(n, e, comments, t, col));
end

save([outputFolder 'result.mat'], 'output');

end
function [isIntersection,nrExits] = detectChiasmataSub(i,nodes,edges,p)
    if mod(i, 100) ==0
        disp(i);
    end
    [thisNodes, thisEdges, thisProb] = connectEM.detectChiasmataPruneToSphere(nodes, edges, ones(size(edges,1),1), p, i);
    C = Graph.findConnectedComponents(thisEdges);
    if length(C) > 3 && sum(cellfun(@(idx)max(pdist2(thisNodes(idx, :), nodes(i,:))) > 4000, C))>3
        isIntersection = true;
        nrExits = length(C);
    else
        isIntersection = false;
        nrExits = 0;
    end
end
