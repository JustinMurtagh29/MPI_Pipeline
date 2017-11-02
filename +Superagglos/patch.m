function newAgglos = patch( ...
        agglo, nodesToAdd, nodesToDelete, ...
        edgesToAdd, edgesToDelete, extraNodeData)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % find nodes to keep
    tempNodeCount = size(agglo.nodes, 1) + size(nodesToAdd, 1);
    nodesToKeep = setdiff(1:tempNodeCount, nodesToDelete);
    
    % update extra node data
    extraNodeData.nodes = nodesToAdd;
    extraNodeFields = fieldnames(extraNodeData);
    
    out = struct;
    for curIdx = 1:numel(extraNodeFields)
        curField = extraNodeFields{curIdx};
        
        curData = agglo.(curField);
        curDataAdd = extraNodeData.(curField);
        
        % sanity check
        assert(size(curData, 1) == size(agglo.nodes, 1));
        assert(size(curDataAdd, 1) == size(nodesToAdd, 1));
        
        % concatenate data
        curData = cat(1, curData, curDataAdd);
        curData = curData(nodesToKeep, :);
        
        out.(curField) = curData;
    end
    
    % renumber node IDs in edges
    edges = cat(1, agglo.edges, edgesToAdd);
    edges(edgesToDelete, :) = [];
    
   [~, edges] = ismember(edges, nodesToKeep);
    assert(all(edges(:)));
    
    % split into connected components
    comps = Graph.findConnectedComponents(edges);
    comps = cat(1, comps, num2cell(reshape(setdiff( ...
        size(out.nodes, 1), cell2mat(comps)), [], 1)));
    
    nodeLUT = Agglo.buildLUT(size(out.nodes, 1), comps);
    edgeComps = nodeLUT(edges(:, 1));
    
    aggloCount = numel(comps);
    newAgglos = struct;
    
    for curIdx = 1:aggloCount
        curNodeIds = comps{curIdx};
        
        % collect node data
        for curFieldIdx = 1:numel(extraNodeFields)
            curField = extraNodeFields{curFieldIdx};
            
            curData = out.(curField);
            curData = curData(curNodeIds, :);
            
            newAgglos(curIdx).(curField) = curData;
        end
        
        % collect edge data
        curEdges = edges(edgeComps == curIdx, :);
        [~, curEdges] = ismember(curEdges, curNodeIds);
        assert(all(curEdges(:)));
        
        newAgglos(curIdx).edges = curEdges;
    end
    
    newAgglos = reshape(newAgglos, [], 1);
end