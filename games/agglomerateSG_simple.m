function [agglo, agglomerated] = agglomerateSG_simple(graph, agglo, t_prob)
    % Agglomerates supervoxel given the graph, always highest probability until under t_prob 

    agglomerated = 0;
    while true
        % Determine neighbours of all agglomerated objects
        theseNeighbours = cat(2,graph.neighbours{agglo});
        % Determine probability of all neighbours
        theseProbabilities = cat(1,graph.neighProb{agglo});
        % Remove self edges (within component)
        selfEdgeIdx = ismember(theseNeighbours, agglo);                      
        theseNeighbours(selfEdgeIdx) = [];
        theseProbabilities(selfEdgeIdx) = [];
        % Determine maxium probability and ID to agglomerated segments
        [maxProb, maxProbIdx] = max(theseProbabilities);
        % ID about to be agglomerated
        maxProbId = theseNeighbours(maxProbIdx);
        % Stop if below probability threshold or over 1000 segments agglomerated
        if maxProb < t_prob || agglomerated >= 5000
            break;
        else
            agglomerated = agglomerated + 1;
            % Add id to agglomeration
            agglo(end+1) = maxProbId;
        end
    end
end

