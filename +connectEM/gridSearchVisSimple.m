% use Manuel's grid search routine with simpler version
% not using heuristics, typeEM information
% visulize the results from the grid search
info = Util.runInfo();

folder = fullfile(p.saveFolder, '20190404T103012gridSearchAgglomeration');
outputFolder = fullfile(folder,'plots');
if ~exist(outputFolder,'dir')
    mkdir(outputFolder)
end
Util.log('Loading results from grid search')
files = dir(fullfile(folder, '*.mat'));
for i=1:length(files)
    result(i) = load(fullfile(folder, files(i).name), ...
        'agglosSize');
end

Util.log('Compare size of mega merger across parameters')
borderSize = [25 50 100 200 300 500]; % border below this size are excluded from dendrite subgraph during inital agglomeration
segmentSize = [25 50 100 300 500 700 1000]; % segment below ...
probThreshold = [0.94 0.95 0.96 0.97 0.98 0.99]; % threshold on neurite continuity probability for CC
%sizeThreshold = 100; % threshold on final agglomerate size in voxels (single segments not collected!)
[a,b,c] = ndgrid(borderSize, segmentSize, probThreshold);
nodesPC = cat(2, a(:), b(:), c(:));
nodesSize = arrayfun(@(x) x.agglosSize(1), result)';
% sanity check
assert(size(nodesSize,1)==size(nodesPC,1))
% sort everything by aggloBig size
[nodesSize, idxSort] = sort(nodesSize,'ascend');
nodesPC = nodesPC(idxSort,:);
% between range 1 and 100 for plotting
nodesSizeMapped = 100*nodesSize./max(nodesSize);
% default colors are sorted in ascending order
nodesColor = jet(size(nodesPC,1));

Util.log('Plotting node Sizes:')
figure; hold on;
plot(nodesSize,'o','MarkerFaceColor','b')
%{out = {};
for i=1:size(nodesPC,1)
    out{i}= sprintf('%d %d %.02f',nodesPC(i,:));
end
allOneString = sprintf(['\''%s\'', '], out{:});
set(gca,'xticklabels',{allOneString(2:end-3)});
xtickangle(90) %}
outfile = fullfile(outputFolder,'aggloSizesDistribution.fig');
export_fig(outfile,'-q101', '-nocrop', '-transparent');


Util.log('Now plotting:')
f = figure();
f.Position = [1 1 21 29.7];
hold on;
a = gca;
for i=1:size(nodesPC,1)
    objectPC = ...
            scatter3(nodesPC(i,1),nodesPC(i,2),nodesPC(i,3), ...
            nodesSizeMapped(i),'filled', 'MarkerEdgeColor', nodesColor(i,:), ...
            'MarkerFaceColor', nodesColor(i,:));
end
colormap('jet')
colorbar
xlabel('bordeSize');
ylabel('segmentSize');
zlabel('probThreshold');
set(gca,'xtick',borderSize)
set(gca,'ytick',segmentSize)
set(gca,'ztick',probThreshold)
grid on
a.XLim = [borderSize(1), borderSize(end)];
a.YLim = [segmentSize(1), segmentSize(end)];
a.ZLim = [probThreshold(1), probThreshold(end)];
outfile = fullfile(outputFolder,'megeMergerSize.fig');
export_fig(outfile,'-q101', '-nocrop', '-transparent');

Util.save(fullfile(outputFolder,'aggloSizes.mat'),result, nodesPC, nodesSize);

