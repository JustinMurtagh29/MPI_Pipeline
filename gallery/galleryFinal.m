function galleryFinal()

skelPath = '/zdata/manuel/sync/fromLap/ekSkel/full/gcl/';

files = dir([skelPath '*.nml']);

for i=1:length(files)
    functionH{i} = @galleryNew;
	inputCell{i} = {skelPath, files(i).name, '/zdata/manuel/sync/wholeCell/retina/gcl_forPetersCorrection/'};
end

startCPU(functionH, inputCell, 'gallery retina');

end

