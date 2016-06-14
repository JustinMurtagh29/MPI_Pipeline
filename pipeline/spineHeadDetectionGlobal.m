function spineHeadDetectionGlobal(p)

me = mfilename;                                            % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));        % where am I located

%Load classifier trained on Kevin's L4 dataset
load([mydir '../spineDetection/shdEnsemble.mat']);

cubeIndices = 1:numel(p.local);
for i = 1:length(cubeIndices)
	cubeNo = cubeIndices(i);
	inputCell{i}={p,cubeNo,Ensemble};
end

functionH = @spineDetectionLocal;

group = ceil(0.10*numel(p.local)); %Group tasks since they are too small to be separately computed on different CPUs

startCPU(functionH,inputCell,'spineHeadDetection',12,group);

end
