function galleryFinal()

skelPath = '/zdata/manuel/sync/fromLap/ekSkel/20150804/all/';
outputDir = '/zdata/manuel/sync/wholeCell/cortex/20150804/gallery/'; 

files = dir([skelPath '*.nml']);
for i=1:length(files)
    functionH{i} = @galleryNew;
	inputCell{i} = {skelPath, files(i).name, outputDir};
end
startCPU(functionH, inputCell, 'gallery retina');

end

