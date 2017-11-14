function output = detectChiasmata(p, nodesV, edges, visualize, outputFolder )
% Detect chiasmata in skeletons based on marching sphere algorithm
% Nodes should be in voxel, scaled here

if ~isfield(p, 'minNrChiasmaExits') ...
        || isempty(p.minNrChiasmaExits)
    % for backward compatibility
    p.minNrChiasmaExits = 4;
end

% Create output folder if it does not exist
if ~isempty(outputFolder) && ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Scale to nm
% Make sure edges are unique
nodes = bsxfun(@times, nodesV, p.voxelSize);
edges = unique(edges, 'rows');

% for each node ("marching sphere" approach to merger detection)
nrExits = zeros(size(nodes, 1), 1);
showProgress = @(n) fprintf('%d of %d nodes done\n', n, size(nodes, 1));

for i = 1:size(nodes, 1)
    nrExits(i) = connectEM.detectChiasmataNodes(p, nodes, edges, i);
    if ~mod(i, 10000); showProgress(i); end
end

% Mark intersections
isIntersection = (nrExits >= p.minNrChiasmaExits);

% Find CC of detected intersections according to graph
[output, queryIdx] = ...
    connectEM.detectChiasmataPostSingleNodeLabel( ...
        p, nodes, nodesV, edges, isIntersection, nrExits);

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
        if any(queryIdx{idx}==-1)
            save([outputFolder 'ohoh'], idx);
            continue;
        end
        n{idx+1} = cat(1, nodesV(output.ccCenterIdx(idx),:), nodesV(output.queryIdx{idx},:));
        e{idx+1} = cat(2, ones(size(n{idx+1},1)-1,1), (2:size(n{idx+1},1))');
        comments{idx+1} = cell(size(n{idx+1},1), 1);
        comments{idx+1}{1} = 'center node of intersection, edges inidcating queries';
        t{idx+1} = ['intersection ' num2str(idx)];
        col{idx+1} = [1 0 0 1];
    end
    writeNml([outputFolder 'skelQueries.nml'], writeSkeletonFromNodesAndEdges(n, e, comments, t, col));
end

if ~isempty(outputFolder)
    fprintf('Writing result... ');
    Util.saveStruct(fullfile(outputFolder, 'result.mat'), output);
    fprintf('done!\n');
end

end
