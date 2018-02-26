% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
wcFile = fullfile(rootDir, 'aggloState', 'wholeCells_08.mat');
outputDir = '/home/amotta/Desktop/whole-cells';

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

wcs = load(wcFile);
wcs = wcs.wholeCells;

%%
mkdir(outputDir);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', mfilename, info.git_repos{1}.hash));

for curIdx = 1:numel(wcs)
    curWc = wcs(curIdx);
    
    curSkel = skel.addTree( ...
        sprintf('Whole cell #%d', curIdx), ...
        curWc.nodes(:, 1:3), curWc.edges);
    curSkel = curSkel.addBranchpoint(find(curWc.axon)); %#ok
    
    curSkelFile = sprintf('whole-cell-%d.nml', curIdx);
    curSkel.write(fullfile(outputDir, curSkelFile));
end
