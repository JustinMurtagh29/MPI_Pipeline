% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
aisFile = fullfile(rootDir, 'aggloState', 'dendrites_17_a.mat');
outDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Load parameter
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% Export AIS super-agglomerates
ais = load(aisFile);
ais = ais.dendrites(ais.indAIS);

% Randomize order
rng(0);
randIds = randperm(numel(ais));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(ais)
    curId = randIds(curIdx);
    curAis = ais(curId);
    
    skelName = sprintf('%d. Axon initial segment %d', curIdx, curId);
    skel = skel.addTree(skelName, curAis.nodes(:, 1:3), curAis.edges);
end

skelFile = fullfile(outDir, 'axon-initial-segments.nml');
skel.write(skelFile);