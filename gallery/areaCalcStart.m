function areaCalcStart()

skelPath = '/gaba/u/mberning/data/retina/retinaN2skeletons/gclNice/';
outputDir = '/gaba/u/mberning/data/retina/retinaN2skeletons/results/surfaceAreaCalculation/'; 

files = dir([skelPath '*.nml']);
for i=1:length(files)
    functionH{i} = @calculateSurfaceInInnerCube;
    inputCell{i} = {skelPath, files(i).name, [outputDir strrep(files(i).name, '.nml', '.mat')]};
end
startCPU(functionH, inputCell, 'surface area calculation');

end

