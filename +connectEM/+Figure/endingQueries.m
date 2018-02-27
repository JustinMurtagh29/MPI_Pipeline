% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_06_c.mat');
outputDir = '/home/amotta/Desktop/ending-queries';

info = Util.runInfo();

%% loading data
axons = load(axonFile);

param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

%% write out examples to webKNOSSOS
bigAxons = axons.indBigAxons;
bigAxons = axons.axons(bigAxons);

axonIds = num2cell(1:numel(bigAxons));
[bigAxons.id] = deal(axonIds{:});

% shuffle axons
rng(0);
sortIds = randperm(numel(bigAxons));
bigAxons = bigAxons(sortIds);

% get rid of hugely merged axon
nodeCount = arrayfun(@(a) size(a.nodes, 1), bigAxons);
bigAxons(nodeCount > 2E3) = [];

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);

skelDesc = sprintf('%s (%s)', info.filename, info.git_repos{1}.hash);
skel = skel.setDescription(skelDesc);

for curIdx = 1:100
    curAxon = bigAxons(curIdx);
    
    curSkel = skel;
    curSkel = curSkel.addTree( ...
        sprintf('Axon %d', curAxon.id), ...
        curAxon.nodes(:, 1:3), curAxon.edges);
    curSkel = curSkel.addBranchpoint( ...
        find(isnan(curAxon.nodes(:, 4)))); %#ok
    
    curSkelFile = sprintf('%02d_axon.nml', curIdx);
    curSkel.write(fullfile(outputDir, curSkelFile));
end

%% selection of axons that illustrate ending queries particularly well
%{
skelDesc = sprintf('%s (%s)', mfilename, info.git_repos{1}.hash);
axonIds = [2369, 4859, 9912, 6285];

bigAxons = axons.indBigAxons;
bigAxons = axons.axons(bigAxons);
bigAxons = bigAxons(axonIds(:));

skels = skeleton();
skels = skels.setDescription(skelDesc);
skels = Skeleton.setParams4Pipeline(skels, param);
skels = Superagglos.buildAggloAndFlightSkels(bigAxons, skels);

for curIdx = 1:numel(skels)
    curSkelFile = sprintf('%d_axon-%d.nml', curIdx, axonIds(curIdx));
    skels(curIdx).write(fullfile(outputDir, curSkelFile));
end
%}