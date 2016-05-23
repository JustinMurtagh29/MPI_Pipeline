function job = saveGlobalSvgDataOnCluster(p)

functionH = @Seg.Global.saveGlobalSvgData;
inputCell = [p];
job = startCPU(functionH,inputCell,'saveGlobalSvgData');

end
