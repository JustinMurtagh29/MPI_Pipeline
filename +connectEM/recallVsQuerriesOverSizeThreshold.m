load '/home/mberning/Desktop/coverage.mat'
queries = load('/home/mberning/Desktop/flightQueries/allQueries.mat');


percentageCovered = arrayfun(@(x)mean(cellfun(@(y)y(1)./y(2), x.axon1.recall_col)), y);
queriesGenerated = arrayfun(@(x)sum(cellfun(@(z)sum(~cat(1,queries.q.exclude{z})), x.axon1.foundAgglomerates_col)), y);

plot(percentageCovered, queriesGenerated);
xlim([0 1]);
ylim([0 360]);
title('Varying size threshold');
xlabel('Node Recall (mean over 10 GT axons)');
ylabel('Number of queries');
