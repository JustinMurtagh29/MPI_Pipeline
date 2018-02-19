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

%%
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

%% write out examples to webKNOSSOS
skelDesc = sprintf('%s (%s)', mfilename, info.git_repos{1}.hash);

for curIdx = 1:50
    curAxon = bigAxons(curIdx);
    
    curSkel = skeleton();
    curSkel = curSkel.addTree( ...
        sprintf('Axon %d', curAxon.id), ...
        curAxon.nodes(:, 1:3), curAxon.edges);
    curSkel = curSkel.addBranchpoint( ...
        find(isnan(curAxon.nodes(:, 4)))); %#ok
    
    curSkel = curSkel.setDescription(skelDesc);
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    
    curSkelFile = sprintf('%02d_axon.nml', curIdx);
    curSkel.write(fullfile(outputDir, curSkelFile));
end