function job = collectGraphStructOnCluster(p)

functionH = @collectGlobalGraphStruct;
inputCell{1} = {p};

job = startCPU(functionH,inputCell,'globalGraphStruct');


end
