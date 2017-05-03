directions.edges = [];
directions.scores = [];
directions.latent = sparse(1, length(segmentMeta.voxelCount));
for idx = 1 : 200
    idx
    directionsPre = load(['/gaba/scratch/kboerg/cluster_directionality_run/', num2str(slice, '%.4u') '.mat']);
    directions.latent(find(directionsPre.y.latent)) = directionsPre.y.latent(find(directionsPre.y.latent));
    directions.edges = [directions.edges; directionsPre.y.edges];
    directions.scores = [directions.scores; directionsPre.y.scores];
end
aggloLookup = [cell2mat(gridAgglo_05{564}.axonsFinal), repelem(1 : length(gridAgglo_05{564}.axonsFinal), cellfun(@length, gridAgglo_05{564}.axonsFinal))'];
[~, ourAgglos] = ismember(find(directions.latent), aggloLookup(:, 1));
directions.agglomerationSize = cellfun(@(x)sum(segmentMeta.voxelCount(gridAgglo_05{564}.axonsFinal(x))));
[~, edgeposition] = ismember(sort(directions.edges, 2), graph.edges, 'rows');
directions.edgeposition = edgeposition;
save([p.saveFolder 'directions.mat'], 'directions')
