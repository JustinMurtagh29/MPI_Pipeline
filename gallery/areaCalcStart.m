function areaCalcStart()

skelPath = '/zdata/manuel/sync/wholeCell/20150505skeletonUpdate/gcl/';

files = dir([skelPath '*.nml']);

for i=1:length(files)
    functionH{i} = @calculateSurfaceInInnerCube;
	inputCell{i} = {skelPath, files(i).name, ['/zdata/manuel/sync/wholeCell/20150505/surfaceAreaCalculation/' strrep(files(i).name, '.nml', '.mat')]};
end

startCPU(functionH, inputCell, 'surface area calculation');

end

