function job = collectGraphStructOnCluster(p)

functionH = @collectGlobalGraphStruct;
inputCell = {p};

job = startCPU(functionH,inputCell,'globalGraphStruct');


end
