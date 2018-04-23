function synAreas = calcSynapseAreas(param, graphT, synT)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxSegId = Seg.Global.getMaxSegId(param);
    
    % Reduce graph
    synLUT = false(maxSegId, 1);
    synLUT(cell2mat(synT.preSegIds)) = true;
    synLUT(cell2mat(synT.postSegIds)) = true;
    
    graphT = graphT(any(synLUT(graphT.edges), 2), :);
    
    % Calculate synapse areas
    synAreas = nan(size(synT, 1), 1);
    for curIdx = 1:size(synT, 1)
        curPreSegIds = synT.preSegIds{curIdx};
        curPostSegIds = synT.postSegIds{curIdx};
        
        % Sanity check
        assert(isempty(intersect( ...
            curPreSegIds, curPostSegIds)));
        
        curMask = ...
            any(ismember(graphT.edges, curPreSegIds), 2) ....
          & any(ismember(graphT.edges, curPostSegIds), 2);
        synAreas(curIdx) = sum(graphT.borderArea(curMask));
    end
end