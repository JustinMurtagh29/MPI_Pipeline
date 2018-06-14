% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
wcFile = fullfile(rootDir, 'aggloState', 'wholeCells_GTAxon_08_v4.mat');
outputDir = '/home/amotta/Desktop/whole-cell-proofreading';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

wcAgglos = load(wcFile);
wcAgglos = wcAgglos.wholeCells;
wcAgglos = SuperAgglo.clean(wcAgglos);

%% 
wcCount = numel(wcAgglos);
numDigits = ceil(log10(1 + wcCount));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:wcCount
    curAgglo = wcAgglos(curIdx);
    curSkel = Superagglos.toSkel(curAgglo, skel);
    
    curSkel.names{1} = sprintf('Whole cell %d', curIdx);
    curSkel = curSkel.addBranchpoint(find(curAgglo.axon)); %#ok
    
    curSkelName = sprintf( ...
        '%0*d_whole-cell-%d.nml', ...
        numDigits, curIdx, curIdx);
    curSkel.write(fullfile(outputDir, curSkelName));
end
