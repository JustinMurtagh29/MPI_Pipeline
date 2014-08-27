function galleryFinal()

skelPath1 = '/zdata/manuel/sync/fromLap/ekSkel/_final/bpc/';
skelPath2 = '/zdata/manuel/sync/fromLap/ekSkel/_final/gcl/';

files1 = dir([skelPath1 '*.nml']);
files2 = dir([skelPath2 '*.nml']);

%matlabpool 12;
for i=1:length(files1)
	galleryNew(skelPath1, files1(i).name, '/zdata/manuel/sync/wholeCell/retina/_final/');
end
for i=1:length(files2)
	galleryNew(skelPath2, files2(i).name, '/zdata/manuel/sync/wholeCell/retina/_final/');
end
%matlabpool close;

end

