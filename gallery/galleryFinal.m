function galleryFinal()

skelPath = '/gaba/u/mberning/data/retina/retinaN2skeletons/allNice/';
outputDir = '/gaba/u/mberning/data/retina/retinaN2skeletons/results/gallery/';

files = dir([skelPath '*.nml']);
for i=1:length(files)
	inputCell{i} = {skelPath, files(i).name, outputDir};
end

functionH = @galleryNew;
startCPU(functionH, inputCell, 'gallery retina');

end

