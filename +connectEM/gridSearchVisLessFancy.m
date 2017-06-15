for idx = 1: 750
    numagglos(idx) = sum(cellfun(@length, result(idx).metrics.axon1.foundAgglomerates_col));
    nummergers(idx) = sum(result(idx).metrics.axon1.mergers_col);
    division=@(x)x(1)/x(2);
    recall(idx) = division(sum(cell2mat(result(idx).metrics.axon1.recall_col')));
    percolators(idx) = result(idx).metrics.axonPercolators(1);
end
subplot(1,2,1)
scatter(numagglos, nummergers)
xlabel('sum agglomerations')
ylabel('sum mergers');
subplot(1,2,2)
scatter(recall, percolators)
xlabel('recall');
ylabel('size of biggest percolator');
figure;
goodones = find(numagglos<312 & nummergers<9);


labels = 'xo'
for idx = 1 : 14
    subplot(1,2,1)
    hold on
scatter(numagglos(goodones(idx))+rand(1)*0.4-0.2, nummergers(goodones(idx))+rand(1)*0.4-0.2,labels(mod(idx,2)+1))
xlabel('sum agglomerations')
ylabel('sum mergers');
title('for 1.7mm total of the 10 axons');
subplot(1,2,2)
hold on
scatter(recall(goodones(idx))+rand(1)*0.004-0.002, percolators(goodones(idx))+rand(1)*1E7-0.5E7,labels(mod(idx,2)+1))

xlabel('recall');
ylabel('size of biggest percolator');
end
legend(strseq('run', goodones))
