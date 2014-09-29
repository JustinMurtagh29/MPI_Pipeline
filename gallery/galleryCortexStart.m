function galleryCortexStart(p)

skelPath = '/zdata/manuel/sync/fromLap/07x2skeletons/wholeCell/';

files = dir([skelPath '*.nml']);

for i=1:length(files)
	galleryCortex(p, skelPath, files(i).name, '/zdata/manuel/sync/wholeCell/retina/_final/');
end

end

