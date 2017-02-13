load('/home/mberning/Desktop/20170203T104102_agglomeration/initialAgglo.mat');

sizePartition = cellfun(@(x)sort(x, 'descend'), sizePartition, 'uni', 0);

% Plot
figure;
hold on;
for i=1:numel(sizePartition)
    plot(sizePartition{i});
end
legend(strseq('Threshold: ', 99:-1:96));
set(gca, 'YScale', 'log');
xlim([0 1e6]);
title('Number elements in equivalence class (sorted)');

x = 99:-1:96;
nrSegmentsAgglomerated = cellfun(@sum, sizePartition)./single(maxSegId); 
volumeAgglomerated = cellfun(@(x)sum(voxelCount(cat(1,x{:}))), initialPartition)./sum(voxelCount);
largetSingleComponent = cellfun(@max, sizePartition);
meanSingleComponent = cellfun(@mean, sizePartition);
figure;
hold on;
subplot(2,2,1);
plot(x, nrSegmentsAgglomerated);
title('%segmentsAgglomerated');
subplot(2,2,2);
plot(x, volumeAgglomerated);
title('%volumeAgglomerated');
subplot(2,2,3);
plot(x, largetSingleComponent);
title('%largestSingleComponentSize');
subplot(2,2,4);
plot(x, meanSingleComponent);
title('%meanSingleComponentSize');



plot(volumeAgglomerated)

plot(largetSingleComponent);

legend('%segmentsAgglomerated', '%volumeAgglomerated', '%sizeLargestSingleComponent', '%segmentsAgglomerated');
set(gca, 'YScale', 'log');
xlim([0 1e6]);
title('Number elements in equivalence class (sorted)');
