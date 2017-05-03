folder = '/gaba/scratch/mberning/aggloGridSearch/';
files = dir([folder 'search06_*.mat']);
for i=1:length(files)
    result(i) = load([folder files(i).name], 'metrics', ...
        'borderSizeDendrites', 'segmentSizeDendrites', 'dendriteProbThreshold', 'probThresholdDendrite', ...
        'borderSizeAxons', 'segmentSizeAxons', 'axonProbThreshold', 'probThresholdAxon');
end
clear i;
outputFolder = '/home/mberning/Desktop/plots4/';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Dendrite viusalization
clearvars -except result files folder outputFolder;

valuesToPlot1 = arrayfun(@(x)cellfun(@(y)y(1)./y(2), x.metrics.dendrite1.recall_col), result, 'uni', 0);
valuesToPlot2 = arrayfun(@(x)x.metrics.dendrite1.mergers_col./x.metrics.dendrite1.length_col*1e6, result, 'uni', 0);
valuesToPlot3 = arrayfun(@(x)(cellfun(@numel, x.metrics.dendrite1.foundAgglomerates_col)-1)./x.metrics.dendrite1.length_col*1e6, result, 'uni', 0);
valuesToPlot4 = arrayfun(@(x)x.metrics.dendritePercolators(1:8), result, 'uni', 0);
valuesToPlot1All = cat(2,valuesToPlot1{:})';
valuesToPlot2All = cat(2,valuesToPlot2{:})';
valuesToPlot3All = cat(2,valuesToPlot3{:})';
valuesToPlot4All = cat(1,valuesToPlot4{:});

groupingSet{1} = arrayfun(@(x)repmat(x.probThresholdDendrite, size(valuesToPlot1{1},2), 1), result, 'uni', 0);
groupingSet{2} = arrayfun(@(x)repmat(x.dendriteProbThreshold, size(valuesToPlot1{1},2), 1), result, 'uni', 0);
groupingSet{3} = arrayfun(@(x)repmat(x.borderSizeDendrites, size(valuesToPlot1{1},2), 1), result, 'uni', 0);
groupingSet{4} = arrayfun(@(x)repmat(x.segmentSizeDendrites, size(valuesToPlot1{1},2), 1), result, 'uni', 0);
groupingSetAll = cat(2, cat(1,groupingSet{1}{:}), cat(1,groupingSet{2}{:}), ...
    cat(1,groupingSet{3}{:}), cat(1,groupingSet{4}{:}));

uniqueValues = unique(groupingSetAll(:,3:4), 'rows');
for i=1:size(uniqueValues,1)
    close all;
    idx = all(bsxfun(@eq, groupingSetAll(:,3:4), uniqueValues(i,:)),2);
    figure('Position', [3841 1 1920 999]);
    subplot(2,2,1);
    boxplot(valuesToPlot1All(idx,:), groupingSetAll(idx,1:2), 'PlotStyle', 'compact');
    set(gca,'FontSize',10,'XTickLabelRotation',90)
    ylabel('Node recall');
    ylim([0.7 1]);
    title(['Dendrites: ' num2str(uniqueValues(i,1)) ' border, ' num2str(uniqueValues(i,2)) ' segment']);
    subplot(2,2,2);
    boxplot(valuesToPlot2All(idx,:), groupingSetAll(idx,1:2), 'PlotStyle', 'compact');
    set(gca,'FontSize',10,'XTickLabelRotation',90)
    ylabel('# merger per mm');
    ylim([0 10]);
    title(['Dendrites: ' num2str(uniqueValues(i,1)) ' border, ' num2str(uniqueValues(i,2)) ' segment']);
    subplot(2,2,3);
    boxplot(valuesToPlot3All(idx,:), groupingSetAll(idx,1:2), 'PlotStyle', 'compact');
    set(gca,'FontSize',10,'XTickLabelRotation',90)
    ylabel('# splits per mm');
    ylim([0 200]);
    title(['Dendrites: ' num2str(uniqueValues(i,1)) ' border, ' num2str(uniqueValues(i,2)) ' segment']);
    subplot(2,2,4);
    boxplot(valuesToPlot4All(idx,:), groupingSetAll(idx,1:2), 'PlotStyle', 'compact');
    set(gca,'FontSize',10,'XTickLabelRotation',90)
    ylabel('Voxel size percolators');
    ylim([1e8 3e10]);
    set(gca, 'YScale', 'log');
    title(['Dendrites: ' num2str(uniqueValues(i,1)) ' border, ' num2str(uniqueValues(i,2)) ' segment']);
    img = getframe(gcf);
    imwrite(img.cdata, [outputFolder 'dendrites_' num2str(uniqueValues(i,1), '%.3i') '_' num2str(uniqueValues(i,2), '%.4i') '.png']);
end

%% Axon visualization
clearvars -except result files folder outputFolder;

valuesToPlot1 = arrayfun(@(x)cellfun(@(y)y(1)./y(2), x.metrics.axon2.recall_col), result, 'uni', 0);
valuesToPlot2 = arrayfun(@(x)x.metrics.axon2.mergers_col./x.metrics.axon2.length_col*1e6, result, 'uni', 0);
valuesToPlot3 = arrayfun(@(x)(cellfun(@numel, x.metrics.axon2.foundAgglomerates_col)-1)./x.metrics.axon2.length_col*1e6, result, 'uni', 0);
valuesToPlot4 = arrayfun(@(x)x.metrics.axonPercolators(1:10), result, 'uni', 0);
valuesToPlot1All = cat(2,valuesToPlot1{:})';
valuesToPlot2All = cat(2,valuesToPlot2{:})';
valuesToPlot3All = cat(2,valuesToPlot3{:})';
valuesToPlot4All = cat(1,valuesToPlot4{:});

groupingSet{1} = arrayfun(@(x)repmat(x.probThresholdAxon, size(valuesToPlot1{1},2), 1), result, 'uni', 0);
groupingSet{2} = arrayfun(@(x)repmat(x.axonProbThreshold, size(valuesToPlot1{1},2), 1), result, 'uni', 0);
groupingSet{3} = arrayfun(@(x)repmat(x.borderSizeAxons, size(valuesToPlot1{1},2), 1), result, 'uni', 0);
groupingSet{4} = arrayfun(@(x)repmat(x.segmentSizeAxons, size(valuesToPlot1{1},2), 1), result, 'uni', 0);
groupingSetAll = cat(2, cat(1,groupingSet{1}{:}), cat(1,groupingSet{2}{:}), ...
    cat(1,groupingSet{3}{:}), cat(1,groupingSet{4}{:}));

uniqueValues = unique(groupingSetAll(:,3:4), 'rows');
for i=1:size(uniqueValues,1)
    close all;
    idx = all(bsxfun(@eq, groupingSetAll(:,3:4), uniqueValues(i,:)),2);
    figure('Position', [3841 1 1920 999]);
    subplot(2,2,1);
    boxplot(valuesToPlot1All(idx,:), groupingSetAll(idx,1:2), 'PlotStyle', 'compact');
    set(gca,'FontSize',10,'XTickLabelRotation',90)
    ylabel('Node recall');
    ylim([0 1]);
    title(['Axons: ' num2str(uniqueValues(i,1)) ' border, ' num2str(uniqueValues(i,2)) ' segment']);
    subplot(2,2,2);
    boxplot(valuesToPlot2All(idx,:), groupingSetAll(idx,1:2), 'PlotStyle', 'compact');
    set(gca,'FontSize',10,'XTickLabelRotation',90)
    ylabel('# merger per mm');
    ylim([0 50]);
    title(['Axons: ' num2str(uniqueValues(i,1)) ' border, ' num2str(uniqueValues(i,2)) ' segment']);
    subplot(2,2,3);
    boxplot(valuesToPlot3All(idx,:), groupingSetAll(idx,1:2), 'PlotStyle', 'compact');
    set(gca,'FontSize',10,'XTickLabelRotation',90)
    ylabel('# splits per mm');
    ylim([0 600]);
    title(['Axons: ' num2str(uniqueValues(i,1)) ' border, ' num2str(uniqueValues(i,2)) ' segment']);
    subplot(2,2,4);
    boxplot(valuesToPlot4All(idx,:), groupingSetAll(idx,1:2), 'PlotStyle', 'compact');
    set(gca,'FontSize',10,'XTickLabelRotation',90)
    ylabel('Voxel size percolators');
    ylim([1e6 1e9]);
    set(gca, 'YScale', 'log');
    title(['Axons: ' num2str(uniqueValues(i,1)) ' border, ' num2str(uniqueValues(i,2)) ' segment']);
    img = getframe(gcf);
    imwrite(img.cdata, [outputFolder 'axons_' num2str(uniqueValues(i,1), '%.3i') '_' num2str(uniqueValues(i,2), '%.4i') '.png']);
end
