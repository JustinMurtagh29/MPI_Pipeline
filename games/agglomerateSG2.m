function [agglo, probabilities] = agglomerateSG2(graph, seeds, nrSteps)
    % Agglomerates supervoxel given the graph, always highest probability for nrSteps

    % Initalization
    agglo = mat2cell(seeds, ones(length(seeds),1));
    probabilities = cell(size(seeds));

    for s=1:length(seeds)
        for n=1:nrSteps
            theseNeighbours = [graph.neighbours{agglo{s}}];
            % Determine probability of all neighbours
            theseProbabilities = cat(1,graph.neighProb{agglo{s}});
            % Remove self edges (within component)
            selfEdgeIdx = ismember(theseNeighbours, agglo{s});
            theseNeighbours(selfEdgeIdx) = [];
            theseProbabilities(selfEdgeIdx) = [];
            % Determine maxium probability and ID to agglomerated segments
            [maxProb, maxId] = max(theseProbabilities);
            agglo{s}(n+1) = theseNeighbours(maxId);
            probabilities{s}(n) = full(maxProb);
        end
    end

end

