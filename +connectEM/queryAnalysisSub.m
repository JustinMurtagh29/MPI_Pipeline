function queryAnalysisSub(agglo,edgesGTall,outputFolder,startidx)


agglo = Superagglos.transformAggloNewOldRepr(agglo);

for idx = startidx : 500 : length(agglo)
    usededges{idx} = find(all(ismember(edgesGTall,agglo{idx}),2));
end
save(fullfile(outputFolder,['output' num2str(startidx)]),'usededges');
