function [newSeeds, probabilities] = agglomerateSG2(graph, seeds, nrSteps)
    % Agglomerates supervoxel given the graph, always highest probability for nrSteps

    % Initalization
    newSeeds = seeds;
    probabilities = cell(size(seeds));
    % Create binary form of adjaceny matrix
    % adjBin = graph.adj ~= 0;

    display('Agglomerating supervoxel');
    tic;
    for s=1:length(seeds)
        for n=1:nrSteps
            % Generate logical vector representation of seed
            % thisSeed = sparse(size(graph.adj,1),1);
            % thisSeed(newSeeds{s}) = 1;
            % Multiply with adjaceny matrix
            % theseProbabilities = graph.adj*thisSeed;
            % Same with binary adjaceny matrix for number edges
            % theseNrEdges = adjBin*thisSeed;
            % Normalize by number of edges
            %theseProbabilities = theseProbabilities./max(theseNrEdges,1);
            % Get maximal edge to each segment
            theseProbabilities = max(bsxfun(@times, graph.adj thisSeed));
            % Remove self edges (within component)
            theseProbabililties = theseProbabilities(newSeeds{s});
            % Add maxium probability ID to agglomerated seed
            [maxProb, maxId] = max(theseProbabilities);
            newSeeds{s}(n+1) = maxId;
            probabilities{s}(n+1) = full(maxProb);
        end
        % Display progress
        Util.progressBar(s, length(seeds));
    end

end

