function galleryFinal()

skelPath1 = '/zdata/manuel/sync/fromLap/ekSkel/_final/bpc/';
skelPath2 = '/zdata/manuel/sync/fromLap/ekSkel/_final/gcl/';

%addpath('/nfs/bmo/mberning/20140310backup/zdataNew/NOBACKUP/code/auxiliary/');
%addpath('/nfs/bmo/mberning/20140310backup/zdataNew/NOBACKUP/code/auxiliary/hocMaker/');
%addpath('/nfs/bmo/mberning/20140310backup/zdataNew/NOBACKUP/code/auxiliary/cubes/');

files1 = dir([skelPath1 '*.nml']);
files2 = dir([skelPath2 '*.nml']);

matlabpool 12;
parfor i=1:length(files1)
	galleryNew(skelPath1, files1(i).name, '/zdata/manuel/sync/wholeCell/retina/_final/');
end
parfor i=1:length(files2)
	galleryNew(skelPath2, files2(i).name, '/zdata/manuel/sync/wholeCell/retina/_final/');
end
matlabpool close;

end

