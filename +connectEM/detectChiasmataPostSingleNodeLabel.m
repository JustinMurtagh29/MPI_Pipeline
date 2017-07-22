function output = detectChiasmataPostSingleNodeLabel(edges, isIntersection, nrExits, nodes, p, nodesV)
temp.edges = edges;
cc = findCCaccordingToGraph(temp, find(isIntersection));
[~, centerOfCC] = cellfun(@(x)min(pdist2(bsxfun(@minus, nodes(x,:), mean(nodes(x,:),1)), [0 0 0])), cc);

% Find out where to query for each CC
queryIdx = cell(length(cc),1);
pos = cell(length(cc),1);
dir = cell(length(cc),1);
dist = cell(length(cc),1);
toDel = false(length(cc),1);
for i=1:length(cc)
    
    % Find all nodes outside CC of detected intersections that are neighbours to CC
    edgesIdx = any(ismember(edges, cc{i}),2);
    queryIdx{i} = setdiff(edges(edgesIdx,:), cc{i});
    % Keep edges not connected to CC
    edgesPruned = edges(~edgesIdx, :);
    nodeDegree = histc(edgesPruned(:), 1:size(nodes,1));
    % Find CC of graph with 5 or more elements after removing CC of intersections
    ccAfterPruning = Graph.findConnectedComponents(edgesPruned, false, true);
    ccAfterPruning = ccAfterPruning(cellfun('length', ccAfterPruning) > 4);
    % Keep only one one queryIdx for each of those CC (the one with maximum node degree)
    queryLocation = cellfun(@(x)intersect(x, queryIdx{i}), ccAfterPruning, 'uni', 0);
    queryLocation = queryLocation(~cellfun('isempty', queryLocation));
    [maxVal, maxIdx] = cellfun(@(x)max(nodeDegree(x)), queryLocation, 'uni', 0);
    if any(cellfun(@(x)x < 1, maxVal))
        % As Graph.findCC is used above with 3rd argument set to true
        error('This should not happen');
    end
    queryIdx{i} = cellfun(@(x,y)x(y), queryLocation, maxIdx);
    % Only generate queries for this CC if at least 2 queries
    if length(queryIdx{i}) > 1
        % Save position direction and length of query
        pos{i} = nodesV(queryIdx{i},:);
        dir{i} = bsxfun(@minus, pos{i}, nodesV(cc{i}(centerOfCC(i)),:));
        dist{i} = pdist2(bsxfun(@times, pos{i}, p.voxelSize), bsxfun(@times, nodesV(cc{i}(centerOfCC(i)),:), p.voxelSize), 'chebychev');
    else
        toDel(i) = true;
    end
end
queryIdx(toDel) = [];
pos(toDel) = [];
dir(toDel) = [];
dist(toDel) = [];
cc(toDel) = [];
centerOfCC(toDel) = [];

% Create an output structure
output.nodes = nodesV;
output.edges = edges;
output.prob = [];
output.isIntersection = isIntersection;
output.nrExits = nrExits;
output.ccNodeIdx = cc;
output.ccCenterIdx = cellfun(@(x,y)x(y), cc, mat2cell(centerOfCC, ones(size(centerOfCC,1),1), ones(size(centerOfCC,2),1)));
output.queryIdx = queryIdx;
output.position = pos;
output.direction = dir;
output.distance = dist;
