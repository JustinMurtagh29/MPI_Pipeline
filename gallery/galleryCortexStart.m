function galleryCortexStart(p)

skelPath = '/zdata/manuel/sync/fromLap/07x2skeletons/wholeCell/';

files = dir([skelPath '*.nml']);

for i=1:length(files)
    functionH{i} = @galleryCortex;
    inputCell{i} = {p, skelPath, files(i).name, '/zdata/manuel/sync/wholeCell/cortex/20141019/'};
end

startCPU(functionH, inputCell, 'whole cell cortex new');

end

