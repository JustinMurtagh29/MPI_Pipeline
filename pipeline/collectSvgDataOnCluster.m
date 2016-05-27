function job = collectSvgDataOnCluster(p)

functionH = @Seg.Global.saveGlobalSvgData;
inputCell = {p};
job = startCPU(functionH,inputCell,'saveGlobalSvgData');

end
