function [agglo, agglomerated] = agglomerateSG_simple(graphStruct, agglo, t_prob, excludeIds)
    % Agglomerates supervoxel given the graph, always highest probability until under t_prob 

    agglomerated = 0;
    while true
        % Determine neighbours of all agglomerated objects
        theseNeighbours = cat(2,graphStruct.neighbours{agglo});
        % Determine probability of all neighbours
        theseProbabilities = cat(1,graphStruct.neighProb{agglo});
        % Remove self edges (within component)
        selfEdgeIdx = ismember(theseNeighbours, agglo);                 
        theseNeighbours(selfEdgeIdx) = [];
        theseProbabilities(selfEdgeIdx) = [];
        % Remove excluded ids
        excEdgeIdx = ismember(theseNeighbours, excludeIds);                      
        theseNeighbours(excEdgeIdx) = [];
        theseProbabilities(excEdgeIdx) = [];
        % Determine maxium probability and ID to agglomerated segments
        [maxProb, maxProbIdx] = max(theseProbabilities);
        % ID about to be agglomerated
        maxProbId = theseNeighbours(maxProbIdx);
        % Stop if below probability threshold or over 10000 segments agglomerated
        if maxProb < t_prob || agglomerated >= 10000
            break;
        else
            agglomerated = agglomerated + 1;
            % Add id to agglomeration
            agglo(end+1,1) = maxProbId;
        end
    end
end

