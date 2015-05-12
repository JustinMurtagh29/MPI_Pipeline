function galleryFinal()

skelPath = '/zdata/manuel/sync/wholeCell/20150505skeletonUpdate/bpc/';

files = dir([skelPath '*.nml']);

for i=1:length(files)
    functionH{i} = @galleryNew;
	inputCell{i} = {skelPath, files(i).name, '/zdata/manuel/sync/wholeCell/20150505/gallery/'};
end

startCPU(functionH, inputCell, 'gallery retina BPC');

clear functionH inputCell;

skelPath = '/zdata/manuel/sync/wholeCell/20150505skeletonUpdate/gcl/';

files = dir([skelPath '*.nml']);

for i=1:length(files)
    functionH{i} = @galleryNew;
	inputCell{i} = {skelPath, files(i).name, '/zdata/manuel/sync/wholeCell/20150505/gallery/'};
end

startCPU(functionH, inputCell, 'gallery retina GCL');

end

