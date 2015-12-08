function [newSeeds, probabilities] = agglomerateSG2(graph, seeds, nrSteps)
    % Agglomerates supervoxel given the graph, always highest probability for nrSteps
    % Only works for single segmentation ID seeds right now

    % Initalization
    newSeeds = seeds;
    probabilities = cell(size(seeds));

    display('Agglomerating supervoxel');
    tic;
    for s=1:length(seeds)
        thisSeed = sparse([],[],[],size(graph.adj,1), size(graph.adj,2), (nrSteps+1)*size(graph.adj,1));
        for n=1:nrSteps
            tic;
            % Generate logical vector representation of seed
            if n==1
                thisSeed(seeds{s}(1),:) = 1;
            else
                thisSeed(maxId,:) = 1;
            end
            % Multiply with GP-based adjaceny matrix
            theseProbabilities = max(graph.adj .* thisSeed,[],1);
            % Remove self edges (within component)
            theseProbabilities(newSeeds{s}) = 0;
            % Add maxium probability ID to agglomerated seed
            [maxProb, maxId] = max(theseProbabilities);
            newSeeds{s}(n+1) = maxId;
            probabilities{s}(n+1) = full(maxProb);
            toc;
        end
        % Display progress
        Util.progressBar(s, length(seeds));
    end

end

