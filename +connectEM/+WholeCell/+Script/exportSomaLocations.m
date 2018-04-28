% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
dendFile  = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03_v2.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

dend = load(dendFile);

%% Generate skeleton
[wcIds, points] = setdiff(dend.idxSomata, 0);

% Reduce somata to points
points = dend.dendrites(points);
points = cell2mat(arrayfun( ...
    @(s) round(mean(s.nodes(:, 1:3), 1)), ...
    points, 'UniformOutput', false));

names = arrayfun( ...
    @(id) sprintf('Whole cell %d', id), ...
    wcIds, 'UniformOutput', false);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

skel = skel.addNodesAsTrees(points, names);
