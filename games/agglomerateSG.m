function newSeeds = agglomerateSG(graph, seeds, threshold, steps, restrictionFlag)
    % Agglomerates supervoxel given the graph, seeds and a probability threshold

    display('Generating efficent thresholded graph representation');
    tic;
    edges = graph.edges(graph.prob > threshold,:);
    maxId = max(edges(:));
    adj = sparse([edges(:,1) edges(:,2)], [edges(:,2) edges(:,1)], 1) + speye(maxId, maxId);
    adjS = adj^steps;
    toc;

    display('Agglomerating supervoxel');
    aggloSeeds = cell(size(seeds));
    tic;
    for i=1:length(seeds)
        % Generate logical vector representation of seed
        thisSeed = sparse(maxId,1);
        thisSeed(seeds{i}) = 1;
        % Multiply with matrix to power of steps (precalculated above)
        temp = adjS*thisSeed;
        aggloSeeds{i} = find(temp);
        Util.progressBar(i, length(seeds));
    end

    if restrictionFlag
        display('Merging redunancies');
        allBorderIds = [seeds{:}];
        cc = sparse(maxId,maxId);
        tic;
        for i=1:length(aggloSeeds)
            temp = intersect(aggloSeeds{i}, allBorderIds);
            [p,q] = meshgrid(temp, temp);
            cc(p(:),q(:)) = 1;
            Util.progressBar(i, length(aggloSeeds));
        end
        % Merge according to connected components
        [eqClasses, objLabels] = findConnectedComponents(cc); 
        newSeeds = eqClasses;
    else
        newSeeds = aggloSeeds;
    end

end

function [equivalenceClasses, objectClassLabels] = findConnectedComponents(sparseAdjMatrixOrEdgeList)
    % Slightly changed version from auxiliaryMethods -> +Graph cause we also want single components
    % Merge there sometime

    % Find block diagonal matrix permutation
    [rowPermutation,~,rowBlockBoundaries] = dmperm(sparseAdjMatrixOrEdgeList);

    % Create vector with one at each row boundary start
    maxValue = size(sparseAdjMatrixOrEdgeList,1);
    newLabelStart = zeros(1,maxValue);
    newLabelStart(rowBlockBoundaries(1:end-1)) = 1;

    % Calculate object class labels (equivalence class of each object)
    objectClassLabels = cumsum(newLabelStart);
    objectClassLabels(rowPermutation) = objectClassLabels;

    % Equivalence classes
    sizeBlocks = diff(rowBlockBoundaries);
    equivalenceClasses = mat2cell(rowPermutation', sizeBlocks);

end

