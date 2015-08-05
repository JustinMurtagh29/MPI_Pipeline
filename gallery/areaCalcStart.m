function areaCalcStart()

skelPath = '/zdata/manuel/sync/fromLap/ekSkel/20150804/gclNice/';
outputDir = '/zdata/manuel/sync/wholeCell/cortex/20150804/surfaceCalculation/'; 

files = dir([skelPath '*.nml']);
for i=1:length(files)
    functionH{i} = @calculateSurfaceInInnerCube;
    inputCell{i} = {skelPath, files(i).name, [outputDir strrep(files(i).name, '.nml', '.mat')]};
end
startCPU(functionH, inputCell, 'surface area calculation');

end

