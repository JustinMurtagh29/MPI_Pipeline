function [agglo, probabilities, excluded] = agglomerateSG6(graph, com, seed, excludeIds, t_prob, t_neigh)
    % Agglomerates supervoxel given the graph, always highest probability until first querry 

    % Initalization
    agglo = seed;
    probabilities = [];
    excluded.ids = [];
    excluded.reasons = {};

    while true
        % Determine neighbours of all agglomerated objects
        theseNeighbours = cat(2,graph.neighbours{agglo});
        % Determine probability of all neighbours
        theseProbabilities = cat(1,graph.neighProb{agglo});
        % Remove self edges (within component)
        selfEdgeIdx = ismember(theseNeighbours, agglo);
        % Remove excluded ids
        excludedEdgeIdx = ismember(theseNeighbours, excluded.ids);
        % Actual removal
        theseNeighbours(selfEdgeIdx | excludedEdgeIdx) = [];
        theseProbabilities(selfEdgeIdx | excludedEdgeIdx) = [];
        % Determine maxium probability and ID to agglomerated segments
        [maxProb, maxProbIdx] = max(theseProbabilities);
        % ID about to be agglomerated
        maxProbId = theseNeighbours(maxProbIdx);
        % Heuristics: Querry skeleton if sth sth sth dark side
        if maxProb < t_prob & length(agglo) < 6
            % If probability is below in first five stops, FORCE start
            agglo(end+1) = maxProbId;
            probabilities(end+1) = maxProb;
        elseif maxProb < t_prob
            % If probability is below t_prob, STOP!
            break;
        elseif any(excludeIds == maxProbId)
            % If id about to be added is in list of tracings, DO NOT AGGLOMERATE
            excluded.ids(end+1) = maxProbId;
            excluded.reasons{end+1} = 'exclusion list';
        elseif length(graph.neighbours{maxProbId}) > t_neigh
            % If id about to be added has more than 40 neighbours, DO NOT AGGLOMERATE
            excluded.ids(end+1) = maxProbId;
            excluded.reasons{end+1} = 'neighbour threshold';
        else
            % If none of the above, actually agglomerate
            agglo(end+1) = maxProbId;
            probabilities(end+1) = maxProb;
        end
    end
end

