function directions = agglomerateDirectionalityCollect(tempfolder, segmentMeta, agglo, graph)
    directions.edges = [];
    directions.scores = [];
    directions.mypca = zeros(length(segmentMeta.voxelCount),9);
    directions.latent = zeros(1, length(segmentMeta.voxelCount));
    for idx = 1 : 200
        idx
        directionsPre = load([tempfolder, num2str(idx, '%.4u') '.mat']);
        directions.latent(find(directionsPre.y.latent)) = directionsPre.y.latent(find(directionsPre.y.latent));
        directions.edges = [directions.edges; directionsPre.y.edges];
        directions.scores = [directions.scores; directionsPre.y.scores];
        directions.mypca(find(directionsPre.y.latent), :) = directionsPre.y.mypca(find(directionsPre.y.latent), :);
    end
    axonsFinalAll = [agglo.axonsFinal,
        num2cell(setdiff(find(segmentMeta.axonProb > 0.5), cell2mat(agglo.axonsFinal)))];
    aggloLookup = [cell2mat(axonsFinalAll), repelem(1 : length(axonsFinalAll), cellfun(@length, axonsFinalAll))'];

    dl = find(directions.latent); % get all segments that have been used
    [~, ourAgglos] = ismember(dl, aggloLookup(:, 1)); % for each segment get the agglo id
    agglomerationSizePre = arrayfun(@(x)sum(segmentMeta.voxelCount(axonsFinalAll{aggloLookup(x,2)})), ourAgglos); %for each segment get the size of the agglo
    directions.agglomerationSize = zeros(1, length(segmentMeta.voxelCount)); %initiate the segmentMeta like array to store the agglo size
    directions.agglomerationSize(dl) = agglomerationSizePre;
    [~, edgeposition] = ismember(sort(directions.edges, 2), graph.edges, 'rows');
    directions.edgeposition = edgeposition;
    save([tempfolder 'directions.mat'], 'directions', '-v7.3')
end
function donothing()
end
