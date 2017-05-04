directions.edges = [];
directions.scores = [];
directions.latent = zeros(1, length(segmentMeta.voxelCount));
for idx = 1 : 200
    idx
    directionsPre = load(['/gaba/scratch/kboerg/cluster_directionality_run/', num2str(idx, '%.4u') '.mat']);
    directions.latent(find(directionsPre.y.latent)) = directionsPre.y.latent(find(directionsPre.y.latent));
    directions.edges = [directions.edges; directionsPre.y.edges];
    directions.scores = [directions.scores; directionsPre.y.scores];
end
axonsFinalAll = [gridAgglo_05{564}.axonsFinal,
    num2cell(setdiff(1 : length(segmentMeta.voxelCount), cell2mat(gridAgglo_05{564}.axonsFinal)))'];
aggloLookup = [cell2mat(axonsFinalAll), repelem(1 : length(axonsFinalAll), cellfun(@length, axonsFinalAll))'];

dl = find(directions.latent); % get all segments that have been used
[~, ourAgglos] = ismember(dl, aggloLookup(:, 1)); % for each segment get the agglo id
agglomerationSizePre = arrayfun(@(x)sum(segmentMeta.voxelCount(axonsFinalAll{aggloLookup(x,2)})), ourAgglos); %for each segment get the size of the agglo
directions.agglomerationSize = zeros(1, length(segmentMeta.voxelCount)); %initiate the segmentMeta like array to store the agglo size
directions.agglomerationSize(dl) = agglomerationSizePre;
[~, edgeposition] = ismember(sort(directions.edges, 2), graph.edges, 'rows');
directions.edgeposition = edgeposition;
save([p.saveFolder 'directions.mat'], 'directions', '-v7.3')
