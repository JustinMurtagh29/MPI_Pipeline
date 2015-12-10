function [newSeeds, probabilities, stats] = agglomerateSG2(graph, seeds, nrSteps)
    % Agglomerates supervoxel given the graph, always highest probability for nrSteps

    % Initalization
    newSeeds = seeds;
    probabilities = cell(size(seeds));
    stats = struct();

    display('Agglomerating supervoxel');
    tic;
    for s=1:length(seeds)
        for n=1:nrSteps
            % Multiply with GP-based adjaceny matrix
            theseProbabilities = max(graph.adj(newSeeds{s},:),[],1);
            % Remove self edges (within component)
            theseProbabilities(newSeeds{s}) = 0;
            % Add maxium probability ID to agglomerated seed
            [maxProb, maxId] = max(theseProbabilities);
            newSeeds{s}(n+1) = maxId;
            probabilities{s}(n+1) = full(maxProb);
            % Temporary section to evaluate stop conditions
            stats(s,n).segId = maxId;
            stats(s,n).neighbours = find(graph.adj(maxId,:));
            stats(s,n).prob = graph.adj(maxId, stats(s,n).neighbours);
            if maxId == 6033884 || maxId ==  6035398 || maxId == 5934237
                break;
            end
        end
        % Display progress
        Util.progressBar(s, length(seeds));
    end

end

