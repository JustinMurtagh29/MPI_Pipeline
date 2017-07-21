function output = detectChiasmata(p, nodesV, edges, visualize, outputFolder )
% Detect chiasmata in skeletons based on marching sphere algorithm
% Nodes should be in voxel, scaled here
%save('/gaba/u/kboerg/biggestpre1.mat','-v7.3');
% Create output folder if it does not exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Scale to nm
% Make sure edges are unique
[edges, ~, idxE] = unique(edges, 'rows');

% Some parameter for algorithm
p.sphereRadiusOuter = Inf; % in nm
p.sphereRadiusInner = 1000; % in nm
p.voxelSize = [11.24 11.24 28];
nodes=bsxfun(@times,nodesV,p.voxelSize);
% for each node ("marching sphere" approach to merger detection)
isIntersection = false(size(nodes,1),1);
nrExits = zeros(size(nodes,1),1);
for i=1:size(nodes,1)
    % if mod(i, 100) ==0
    %     disp(i);
    % end
    [thisNodes, thisEdges, thisProb] = connectEM.detectChiasmataPruneToSphere(nodes, edges, ones(size(edges,1),1), p, i);
    C = Graph.findConnectedComponents(thisEdges);
    if length(C) > 3 && sum(cellfun(@(idx)max(pdist2(thisNodes(idx, :), nodes(i,:))) > 4000, C))>3
        isIntersection(i) = true;
        nrExits(i) = length(C);
    end
    
    if visualize && isIntersection(i)
        figure('Position', [3841 1 1920 999]);
        % First subplot visualizing pruning to sphere (Step 1)
        subplot(2,2,1);
        visualizeSingleSphere(thisNodes, thisEdges, thisProb, thisIdx, thisIdxOuter, p);
        colorbar;
        title('Pruned to CC of edges within outer sphere');
        % Second subplot visualizing result of inital clustering
        subplot(2,2,2);
        visualizeClustering(thisNodes, thisEdges, thisProb, thisIdx, thisIdxOuter, p, thisNodes3, bestClustering.idxAll);
        colorbar;
        title('Clusters in sphere hull as detected by cluster visualization');
        % Third subplot
        subplot(2,2,3);
        imagesc(distances);
        axis equal; axis off;
        colorbar;
        title('Cosine distances between nodes in outer sphere');
        subplot(2,2,4);
        imagesc(distances > .1);
        axis equal; axis off;
        colorbar;
        title(['Final result: Detected intersection with ' num2str(nrExits(i)) ' exits']);
        % Save both as fig and png
        saveas(gcf, [outputFolder num2str(i) '.fig']);
        %img = getframe(gcf);
        %imwrite(img.cdata, [outputFolder num2str(i) '.png']);
        close all;
    end
end
%save('/gaba/u/kboerg/biggest1.mat', 'isIntersection');
% Find CC of detected intersections according to graph
output = detectChiasmataPostSingleNodeLabel(edges, isIntersaction, nrExits, nodes, p);

if visualize
    % Write result to skletons for control (detection of intersections)
    comments = cell(size(nodesV,1), 1);
    comments(isIntersection) = arrayfun(@(x)strcat('intersection, ', num2str(x), ' exits'), nrExits(isIntersection), 'uni', 0);
    writeNml([outputFolder 'skelDetected.nml'], writeSkeletonFromNodesAndEdges({nodesV}, {edges}, {comments}, {['axon_' num2str(i, '%.5i')]}, {[0 0 1 1]}));
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
    for i=1:length(queryIdx)
        n{i+1} = cat(1, nodesV(output.ccCenterIdx(i),:), nodesV(output.queryIdx{i},:));
        e{i+1} = cat(2, ones(size(n{i+1},1)-1,1), (2:size(n{i+1},1))');
        comments{i+1} = cell(size(n{i+1},1), 1);
        comments{i+1}{1} = 'center node of intersection, edges inidcating queries';
        t{i+1} = ['intersection ' num2str(i)];
        col{i+1} = [1 0 0 1];
    end
    writeNml([outputFolder 'skelQueries.nml'], writeSkeletonFromNodesAndEdges(n, e, comments, t, col));
end

save([outputFolder 'result.mat'], 'output');

end



function [thisNodes, thisEdges, thisProb, thisIdx, thisIdxOuter] = restrictGraphToCC(thisNodes, thisEdges, thisProb, thisIdx, thisIdxOuter)

% Keep only CC connected to "current" node
thisCC = Graph.findConnectedComponents(thisEdges, false, false);
thisCC = thisCC{cellfun(@(x)any(ismember(x,thisIdx)), thisCC)};
thisNodes = thisNodes(thisCC,:);
thisIdx = find(ismember(thisCC,thisIdx));
thisIdxOuter = intersect(thisIdxOuter, thisCC);
for i=1:length(thisIdxOuter)
    thisIdxOuter(i) = find(ismember(thisCC,thisIdxOuter(i)));
end
tempIdx = all(ismember(thisEdges,thisCC),2);
thisEdges = thisEdges(tempIdx,:);
thisEdges = arrayfun(@(x)find(x == thisCC), thisEdges);
thisProb = thisProb(tempIdx);

end





function visualizeSingleSphere(nodes, edges, prob, currentNodeIdx, outerNodesIdx, p)

rP = p.sphereRadiusOuter;
rC = p.sphereRadiusInner;

% Define some values
nrBins = 10;
bins = linspace(0, 1, nrBins+1);

% RGB colors for indicating edge weight
edgeColors(:,1) = linspace(1,0,nrBins);
edgeColors(:,2) = linspace(0,1,nrBins);
edgeColors(:,3) = zeros(1,10);
% Calculate color for each edge
[row, col] = find(bsxfun(@gt, prob, bins(1:end-1)) & bsxfun(@le, prob, bins(2:end)));
[~, rIdx] = sort(row);
thisProbBinned = col(rIdx);
thisEdgeColors = edgeColors(thisProbBinned,:);

% Plot
hold on;
for i=1:size(edges,1)
    plot3([nodes(edges(i,1),1) nodes(edges(i,2),1)]', ...
        [nodes(edges(i,1),2) nodes(edges(i,2),2)]', ...
        [nodes(edges(i,1),3) nodes(edges(i,2),3)]', 'Color', thisEdgeColors(i,:));
end
plot3(nodes(currentNodeIdx, 1), nodes(currentNodeIdx, 2),nodes(currentNodeIdx, 3), 'xb', 'MarkerSize', 10);
plot3(nodes(outerNodesIdx, 1), nodes(outerNodesIdx, 2),nodes(outerNodesIdx, 3), 'xb', 'MarkerSize', 10);
% Add spheres
[x,y,z] = sphere;
surf(rP*x+nodes(currentNodeIdx,1),rP*y+nodes(currentNodeIdx,2),rP*z+nodes(currentNodeIdx,3), 'EdgeColor', 'none', 'FaceColor', 'b');
surf(rC*x+nodes(currentNodeIdx,1),rC*y+nodes(currentNodeIdx,2),rC*z+nodes(currentNodeIdx,3), 'EdgeColor', 'none', 'FaceColor', 'b');
alpha(.2);
camlight;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view(3);
colormap(edgeColors);
caxis([0 1]);

end

function visualizeClustering(nodes, edges, prob, currentNodeIdx, outerNodeIdx, p, clusterNodes, clusterIdx)

clusterNodes = bsxfun(@plus, clusterNodes, nodes(currentNodeIdx,:));

hold on;
visualizeSingleSphere(nodes, edges, prob, currentNodeIdx, outerNodeIdx, p);
clusters = unique(clusterIdx);
markers = {'or' 'og' 'ob' 'oy'};
for i=1:length(clusters)
    tempIdx = clusters(i) == clusterIdx;
    plot3(clusterNodes(tempIdx,1), clusterNodes(tempIdx,2), clusterNodes(tempIdx,3), markers{i}, 'MarkerSize', 10);
end

end
