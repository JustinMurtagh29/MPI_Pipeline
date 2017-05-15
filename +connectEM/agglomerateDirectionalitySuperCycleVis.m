for idx = 1:2
    cycleAgglo_01{idx} = load(['/gaba/scratch/kboerg/directCycle/cycle' num2str(idx, '%0.3u')])
end
for idx = 1 : 3
    figure
    a=[1, 2, 4, 8, 16, 32, 64, 128, 256, 512];
    [N,edges] = histcounts(cellfun('length', cycleAgglo_01{idx}.axonsFinal),[sqrt(a(1:end-1).*a(2:end)), Inf]);
    bar(log10(N)+0.5);
    ylabel('counts')
    set(gca, 'YTick', 0.5:6.5, 'YTickLabel', {'1','10','100','1k','10k','100k','1M'});
    ylim([0.2,6.7]);
    set(gca, 'XTick', 1:11, 'XTickLabel', {'2^1','2^2','2^3','2^4','2^5','2^6','2^7','2^8','2^9'});
    saveas(gcf, sprintf('/gaba/scratch/kboerg/directCycle/figure%u.pdf', idx));

end
