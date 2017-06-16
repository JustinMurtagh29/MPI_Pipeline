load /home/mberning/Desktop/flightQueries/allQueries.mat;

% Plot cumulative histogram for MH
queriesNotExcluded = cellfun(@(x)sum(~x), q.exclude);
voxelCount = cell2mat(q.voxelCount);
edges = linspace(quantile(voxelCount, 0.0001), quantile(voxelCount, 0.9), 10002);
voxelCountTemp = voxelCount;
voxelCountTemp(voxelCountTemp >= edges(end)) = edges(end)-eps;
voxelCountTemp(voxelCountTemp <= edges(1)) = edges(1)+eps;
queryHist = histcounts(voxelCountTemp, edges);
queryCum = cumsum(queryHist(end:-1:1));
queryEdges = arrayfun(@(x,y)round(mean(cat(1,x,y))), edges(end:-1:2), edges(end-1:-1:1));

figure;
plot(queryCum);
set(gca, 'XTick', 1:500:numel(queryCum));
set(gca, 'XTickLabel', arrayfun(@(x)num2str(x), queryEdges(1:500:numel(queryCum)), 'uni', 0));
set(gca, 'YScale', 'log');
xlim([1 numel(queryCum)]);
xlabel('Agglomerate size [voxel]');
ylabel('# Queries in agglomerates above this size');
