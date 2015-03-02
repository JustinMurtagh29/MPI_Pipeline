function galleryCortexStart(p)

skelPath = '/zdata/manuel/sync/fromLap/07x2skeletons/axonsMHforPaper/';
%skelPath = '/zdata/manuel/sync/fromLap/07x2skeletons/spinyStellateForPaper/';

files = dir([skelPath '*.nml']);

for i=1:length(files)
    functionH{i} = @galleryCortexSingleIsosurfaces;;
    inputCell{i} = {p, skelPath, files(i).name, '/zdata/manuel/sync/wholeCell/cortex/20150301/'};
end

startCPU(functionH, inputCell, 'whole cell cortex new');

end

