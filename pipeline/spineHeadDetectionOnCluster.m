function job =spineHeadDetectionOnCluster(p)

cubeIndices = 1:numel(p.local);
for i = 1:length(cubeIndices)
	cubeNo = cubeIndices(i);
	inputCell{i}={p,cubeNo};
end

functionH = @spineHeadDetectionLocal;

group = ceil(0.10*numel(p.local)); %Group tasks since they are too small to be separately computed on different CPUs

job = startCPU(functionH,inputCell,'spineHeadDetection',12,group);
display('Now calculating spine heads in local cubes');
end
