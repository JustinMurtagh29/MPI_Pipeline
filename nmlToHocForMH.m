addpath(genpath('/home/mberning/code/SegEM/auxiliaryMethods'));
addpath(genpath('/home/mberning/code/manuelCode/visualization'));

folder = '/home/mberning/Desktop/forMH/';

files = dir([folder '*.nml']);
for i=1:length(files)
    skel = parseNml([folder files(i).name]);
    %skel{1}.nodes(:,1:3) = bsxfun(@times, skel{1}.nodes(:,1:3), [12 12 25]);
    convertKnossosNmlToHoc2(skel, [folder files(i).name(1:end-4)], 0, 1, 0, 0, [11.24 11.24 28]);
end