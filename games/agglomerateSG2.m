function [newSeeds, probabilities] = agglomerateSG2(graph, seeds, nrSteps)
    % Agglomerates supervoxel given the graph, always highest probability for nrSteps

    % Initalization
    newSeeds = seeds;
    probabilities = cell(size(seeds));

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
        end
        % Display progress
        Util.progressBar(s, length(seeds));
    end

end

