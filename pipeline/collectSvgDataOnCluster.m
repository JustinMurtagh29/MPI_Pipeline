function job = collectSvgDataOnCluster(p)

functionH = @Seg.Global.saveGlobalSvgData;
inputCell{1} = {p};
job = startCPU(functionH,inputCell,'saveGlobalSvgData',[],[],-100);

end
