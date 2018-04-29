% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
outputDir = '/home/amotta/Desktop/whole-cells';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

dend = load(dendFile);

%% Generate skeletons
[~, aggloIds] = setdiff(dend.idxWholeCells, 0);
agglos = dend.dendrites(aggloIds);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

numDigits = ceil(log10(1 + numel(agglos)));

for curIdx = 1:numel(agglos)
    curId = aggloIds(curIdx);
    curAgglo = agglos(curIdx);
    
    curSkel = connectEM.Tweak.buildNml(curAgglo, curId, skel);
    curSkelName = sprintf('%0*d_agglo-%d.nml', numDigits, curIdx, curId);
    curSkel.write(fullfile(outputDir, curSkelName));
end
