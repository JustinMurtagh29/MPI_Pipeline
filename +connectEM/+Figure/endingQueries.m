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
%{
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
%}

%% selection of axons that illustrate ending queries particularly well
skelDesc = sprintf('%s (%s)', mfilename, info.git_repos{1}.hash);
axonIds = [2369, 4859, 9912, 6285];

bigAxons = axons.indBigAxons;
bigAxons = axons.axons(bigAxons);
bigAxons = bigAxons(axonIds(:));

for curIdx = 1:numel(axonIds)
    curSagglo = bigAxons(curIdx);
    assert(issorted(curSagglo.edges, 2));
    
    % separate flight paths from agglomerates
    curIntraEdges = ismember( ...
        curSagglo.edges, find(isnan(curSagglo.nodes(:, 4))));
    curIntraEdges = curSagglo.edges(sum(curIntraEdges, 2) ~= 1, :);
    
   [curNodeIds, ~, curIntraEdges] = unique(curIntraEdges);
    curIntraEdges = reshape(curIntraEdges, [], 2);
    
    curNodeCount = numel(curNodeIds);
    curNodes = curSagglo.nodes(curNodeIds, :);
    
    % build agglomerates
    curAdj = sparse( ...
        curIntraEdges(:, 2), curIntraEdges(:, 1), ...
        true, curNodeCount, curNodeCount);
   [curCompCount, curLUT] = ...
        graphconncomp(curAdj, 'Directed', false);
    assert(all(curLUT > 0));
    
    % split super-agglomerate
    curSkel = skeleton;
    for curCompIdx = 1:curCompCount
        curCompNodeIds = find(curLUT == curCompIdx);
        curCompNodes = curNodes(curCompNodeIds, 1:3);
        
        % edges
       [~, curCompEdges] = ismember(curIntraEdges, curCompNodeIds);
        curCompEdges(~all(curCompEdges, 2), :) = [];
        
        if any(isnan(curNodes(curCompNodeIds, 4)))
            % flight path
            curCompName = 'Flight';
            curCompColor = [1, 0, 0, 1];
        else
            % agglomerate
            curCompName = 'Agglomerate';
            curCompColor = [0, 0, 1, 1];
        end
        
        % name
        curCompName = sprintf( ...
            '%0*d %s', ceil(log10(curCompCount + 1)), ...
            curCompIdx, curCompName);
        curSkel = curSkel.addTree( ...
            curCompName, curCompNodes, curCompEdges, curCompColor);
    end
    
    curSkel = curSkel.setDescription(skelDesc);
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    
    curSkelFile = sprintf('%d_axon-%d.nml', curIdx, axonIds(curIdx));
    curSkel.write(fullfile(outputDir, curSkelFile));
end