%%
a=y6cycle;
b=yAM;

close all
subplot(1,3,1)
title('b black');
rng(20170804)
hold on
labels = 'xo'
for idx = 1 : 10    
    firstrand(idx) = 0.2*rand(1)-0.1
    plot([b.axon1.length_col(idx)/1000; a.axon1.length_col(idx)/1000]', [cellfun(@length,b.axon1.foundAgglomerates_col(idx))+firstrand(idx); cellfun(@length,a.axon1.foundAgglomerates_col(idx))+0.2*rand(1)-0.1]', ['-' labels(mod(idx,2)+1)]);
end
ylabel('# of agglomerates')
xlabel('length [µm]')
legend(strseq('axon ',1:10))
subplot(1,3,2)
hold on
for idx = 1 : 10
    plot([b.axon1.length_col(idx)/1000; a.axon1.length_col(idx)/1000]', [cellfun(@(x)x(1)/x(2),b.axon1.recall_col(idx)); cellfun(@(x)x(1)/x(2),a.axon1.recall_col(idx))]', ['-' labels(mod(idx,2)+1)]);
end
ylabel('fraction covered')
xlabel('length [µm]')


subplot(1,3,3)
hold on
for idx = 1: 10
    plot([cellfun(@length,b.axon1.foundAgglomerates_col(idx))+firstrand(idx);cellfun(@length,a.axon1.foundAgglomerates_col(idx))+0.2*rand(1)-0.1],[b.axon1.mergers_col(idx)+firstrand(idx);a.axon1.mergers_col(idx)+0.2*rand(1)-0.1], ['-' labels(mod(idx,2)+1)]);
end
ylabel('mergers')
xlabel('# of agglomerates')
subplot(1,3,1)
for idx = 1 : 10
    scatter(b.axon1.length_col(idx)/1000, cellfun(@length,b.axon1.foundAgglomerates_col(idx))+firstrand(idx), ['k' labels(mod(idx,2)+1)]);
end
subplot(1,3,2)
for idx = 1 : 10
    scatter(b.axon1.length_col(idx)/1000, cellfun(@(x)x(1)/x(2),b.axon1.recall_col(idx)), ['k' labels(mod(idx,2)+1)]);
end
subplot(1,3,3)
for idx = 1: 10
    scatter(cellfun(@length,b.axon1.foundAgglomerates_col(idx))+firstrand(idx),b.axon1.mergers_col(idx)+firstrand(idx), ['k' labels(mod(idx,2)+1)]);
end
