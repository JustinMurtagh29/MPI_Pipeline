function out = applyNml(agglo, nml)
    % agglo = applyNml(agglo, nml)
    %   Applies manual changes to the super-agglomerates.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % Prepare node ID translation table
    comments = NML.buildCommentTable(nml);
    
    comments.origId = regexpi( ...
        comments.comment, 'Node\s+(?<id>\d+)', 'names', 'once');
    comments.origId = cellfun(@(n) str2double(n.id), comments.origId);
    
    % Sanity checks
    assert(numel(comments.node) ...
        == numel(unique(comments.node)));
    assert(numel(comments.origId) ...
        == numel(unique(comments.origId)));
    
    newNodeCount = max(cellfun( ...
        @(nodes) max(nodes.id), nml.things.nodes));
    origIdLUT = zeros(newNodeCount, 1);
    origIdLUT(comments.node) = comments.origId;

    %% Build new output agglomerates
    out = struct('nodes', {}, 'edges', {});
    trees = NML.buildTreeTable(nml);

    for curIdx = 1:size(trees, 1)
        curNodes = trees.nodes{curIdx};
        curEdges = trees.edges{curIdx};

        % Nodes
        curAgglo = struct;
        curAgglo.nodes = [ ...
            1 + curNodes.x, ...
            1 + curNodes.y, ...
            1 + curNodes.y, ...
            nan(size(curNodes.id))];

        % Re-use old data
        curOrigNodeIds = curNodes.id;
        curOrigNodeIds = origIdLUT(curOrigNodeIds);

        curMask = curOrigNodeIds ~= 0;
        curAgglo.nodes(curMask, :) = ...
            agglo.nodes(curOrigNodeIds(curMask), :);

        % Edges
        curAgglo.edges = [ ...
            curEdges.source, ...
            curEdges.target];
       [~, curAgglo.edges] = ismember( ...
            curAgglo.edges, curNodes.id);

        out(curIdx) = curAgglo;
    end

    out = reshape(out, [], 1);
end
