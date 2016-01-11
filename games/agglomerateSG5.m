function [agglo, probabilities, mergerList, queryId] = agglomerateSG5(graph, com, seeds, t_prob, t_neigh)
    % Agglomerates supervoxel given the graph, always highest probability until first querry 

    % Initalization
    agglo = mat2cell(seeds, ones(length(seeds),1));
    probabilities = cell(size(agglo));
    querry = false;
    mergerList = [];

    while ~querry
        % Remove merger from agglomerated segments
        thisAggloWithoutMerger = setdiff(agglo{1}, mergerList);
        % Determine neighbours of all agglomerated objects
        theseNeighbours = [graph.neighbours{thisAggloWithoutMerger}];
        % Determine probability of all neighbours
        theseProbabilities = cat(1,graph.neighProb{thisAggloWithoutMerger});
        % Remove self edges (within component)
        selfEdgeIdx = ismember(theseNeighbours, agglo{1});                      
        theseNeighbours(selfEdgeIdx) = [];
        theseProbabilities(selfEdgeIdx) = [];
        % Determine maxium probability and ID to agglomerated segments
        [maxProb, maxProbIdx] = max(theseProbabilities);
        % ID about to be agglomerated
        maxProbId = theseNeighbours(maxProbIdx);
        % Heuristics: Querry skeleton if sth sth sth dark side
        if maxProb < t_prob
            % Determine segment farest way from start com
            startCoord = com(seeds(1),:);
            aggloCoord = com(agglo{1},:);
            distances = sqrt(sum(bsxfun(@times, bsxfun(@minus, aggloCoord, startCoord), [11.24 11.24 28]).^2,2));
            [~,maxDistIdx] = max(distances);
            % Querry skeleton
            queryId = agglo{1}(maxDistIdx);
            break;
            querry = true;
        else
            addedId = maxProbId;
        end
        % Add id to agglomeration
        agglo{1}(end+1) = addedId;
        probabilities{1}(end+1) = maxProb;
        % Next Heuristics: Exclude if 40 neighbours or more
        if length(graph.neighbours{addedId}) >= t_neigh
            % Add last seed ID to ignore list (SegEM merger)
            mergerList(end+1) = addedId;
        end
    end
end

