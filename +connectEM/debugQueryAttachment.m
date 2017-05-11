function debugQueryAttachment(segmentPositions, agglos, ff, outputFolder)

    % Visualize each as one nml (set of trees)
    for i=1:length(ff.segIds)
        % Initalize variables
        nodes = {};
        treeNames = {};
        comments = {};
        % Determine occurence of each segment in node + neighbours (seperate start node)
        [uniqueSegments, neighboursStartNode, nodesExcludedIdx, startNodeIdx] = connectEM.queryAnalysis( ...
            ff.segIds{i}, ff.neighbours{i,:}, ff.nodes{i}, ff.startNode{i}); 
        % Determine which agglomerations overlap with given query 
        overlapWithAgglomerations = cellfun(@(x)intersect(x, uniqueSegments.segIds), agglos, 'uni', 0);
        overlapsFound = find(~cellfun('isempty', overlapWithAgglomerations));
        % For each node in the inital query, add comment if needed (bounding box, start node, or excluded or leftover segments)
        nodes{1} = ff.nodes{i};
        comments{1}(nodesExcludedIdx) = {'out of bounding box'};
        startNodeComment = 'start node';
        for j=1:length(neighboursStartNode.segIds)
            startNodeComment = strcat(startNodeComment, ', ', num2str(neighboursStartNode.segIds(j)), ': ', num2str(neighboursStartNode.occurences(j)));
        end
        comments{1}(startNodeIdx) = {startNodeComment};
        treeNames{1} = 'Query';
        % For each agglomeration found add it to nodes & add statistics in treeName 
        for j=1:length(overlapsFound)
            nodes{1+j} = segmentPositions(agglos{overlapsFound(j)},:);
            tempIdx = ismember(uniqueSegments.segIds, overlapWithAgglomerations{overlapsFound(j)});
            tempVal = sum(uniqueSegments.occurences(tempIdx));
            treeNames{1+j} = ['Agglomerate ' num2str(overlapsFound(j), '%.2i'), ', Node evidence: ' num2str(tempVal, '%.3i')];
        end       
        % Compose & write skeletons for visualization
        nodes = cellfun(@(x)bsxfun(@minus, x, [1 1 1]), nodes, 'uni', 0);
        connectEM.generateSkeletonFromNodes([outputFolder 'query_' num2str(i, '%.3i') '.nml'], nodes, treeNames, comments);
    end

end

