addpath(genpath('/u/kostalr/code/'));

root = 'u/kostalr/data/data/thirdcubing/4/';
prefix = '2016-03-22_ex4_RN_CC_S1BF_L4_mag4';

raw = readKnossosRoi(root, prefix, [160 1445; 192 1860; 32 930]);

% Detect vessels with visualization flag set to false
% Set flag to true for new dataset parameter tuning
vessels = detectVessels(raw, false);
save('/u/kostalr/newtry/vessels/vs.mat','vessels');
% Detect nuclei with visualization flag set to false
% Set flag to true for new dataset parameter tuning
%nuclei = detectNuclei(raw, vessels, false);

% Dilate with sphere to get outside nucleus and capture at least some endothelial segments
[x,y,z] = meshgrid(-3:3,-3:3,-1:1);
se = (x/3).^2 + (y/3).^2 + (z/1).^2 <= 1;
nuclei = imdilate(nuclei, se);
vessels = imdilate(vessels, se);

% Make movie and isosurface visualization of detected vessels
makeSegMovie(vessels, raw, '/u/kostalr/newtry/vessels.avi', 1);



% Make movie and isosurface visualization of detected nuclei
%makeSegMovie(nuclei, raw, '/u/kostalr/newtry/nuclei.avi', 1);



nucleiL = labelmatrix(bwconncomp(nuclei));
rp = regionprops(nucleiL);


writeKnossosRoi('/u/kostalr/newtry/vessels/', 'vessels', [160 192 32], vessels);
%writeKnossosRoi('/u/kostalr/newtry/nuclei/', 'nuclei', [160 192 32], nucleiL);
