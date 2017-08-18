function queryAnalysisSub(agglo,edgesGTall,outputFolder,startidx)

if isfield(agglo,'nodes') % transform new representation into old
    agglo = cellfun(@(x) x(:,4),{agglo.nodes},'uni',0);
end
for idx = startidx : 500 : length(agglo)
    usededges{idx} = find(all(ismember(edgesGTall,agglo{idx}),2));
end
save(fullfile(outputFolder,['output' num2str(startidx)]),'usededges');
