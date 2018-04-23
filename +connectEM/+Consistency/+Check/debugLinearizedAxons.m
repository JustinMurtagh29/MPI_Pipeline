% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
origFile = fullfile(rootDir, 'aggloState', 'axons_18_b.mat');
splitFile = fullfile(rootDir, 'aggloState', 'axons_18_b_linearized.mat');

nmlDir = '/home/amotta/Desktop/linearized';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

orig = load(origFile);
split = load(splitFile);

%% Export examples
splitIds = accumarray(split.parentIds, 1);
splitIds = find(splitIds > 1);

rng(0);
randIds = randperm(numel(splitIds));
randIds = splitIds(randIds);
randIds = randIds(1:25);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

numDigits = ceil(log10(1 + numel(randIds)));

for curIdx = 1:numel(randIds)
    curId = randIds(curIdx);
    
    curOrig = orig.axons(curId);
    curSplit = split.parentIds == curId;
    curSplit = split.axons(curSplit);
    
    curSkel = skel;
    curSkel = Superagglos.toSkel(curOrig, curSkel);
    curSkel = Superagglos.toSkel(curSplit, curSkel);
    
    curNames = sprintf('Axon %d', curId);
    curNames = [curNames; arrayfun( ...
        @(id) sprintf('Part %d', id), ...
        transpose(1:numel(curSplit)), ...
        'UniformOutput', false)]; %#ok
    curSkel.names = curNames;
    
    curSkelName = sprintf('%0*d_axon-%d.nml', numDigits, curIdx, curId);
    curSkel.write(fullfile(nmlDir, curSkelName));
end
