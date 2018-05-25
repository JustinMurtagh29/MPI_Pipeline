% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
wcFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');

outCellIds = [5, 20, 25, 40, 54, 66, 83, 86];
outDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

wcData = load(wcFile);

%% Export super-agglomerates
numDigits = ceil(log10(1 + numel(outCellIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(outCellIds)
    curCellId = outCellIds(curIdx);
    
    curAgglo = wcData.idxWholeCells == curCellId;
    curAgglo = wcData.dendrites(curAgglo);
    
    curSkel = Superagglos.toSkel(curAgglo, skel);
    curSkel.names{end} = sprintf('Whole cell %d', curCellId);
    
    curSkelName = fullfile(outDir, sprintf( ...
        '%0*d_whole-cell-%d.nml', numDigits, curIdx, curCellId));
    curSkel.write(curSkelName);
end
