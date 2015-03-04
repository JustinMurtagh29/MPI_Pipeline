function galleryCortexStart(p)

%skelPath = '/zdata/manuel/sync/fromLap/07x2skeletons/axonsMHforPaper/';
%skelPath = '/zdata/manuel/sync/fromLap/07x2skeletons/spinyStellateForPaper/';
skelPath = '/zdata/manuel/sync/fromLap/07x2skeletons/forDendritesChapter/';
outputPath = '/zdata/manuel/sync/wholeCell/cortex/20150301/';

if ~exist(outputPath)
    mkdir(outputPath);
end

files = dir([skelPath '*.nml']);

for i=1:length(files)
    functionH{i} = @galleryCortexSingleIsosurfaces;;
    inputCell{i} = {p, skelPath, files(i).name, outputPath};
end

startCPU(functionH, inputCell, 'whole cell cortex new');

end

