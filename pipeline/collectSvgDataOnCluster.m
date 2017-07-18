function job = collectSvgDataOnCluster(p,runlocal)

if nargin < 2
	runlocal = 0;
end

if runlocal
	Seg.Global.saveGlobalSvgData(p,[],[],1);
	job = [];
else
functionH = @Seg.Global.saveGlobalSvgData;
inputCell{1} = {p};
job = startCPU(functionH,inputCell,'saveGlobalSvgData',[],[],-100);
end

end
