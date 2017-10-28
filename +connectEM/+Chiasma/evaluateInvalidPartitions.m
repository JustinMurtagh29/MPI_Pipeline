% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_08_a.mat');
outputDir = '/home/amotta/Desktop/unsolved-chiasmata';

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

data = load(axonFile);

%% find unsolved and invalid chiasmata chiasmata
chiasmata = table;
chiasmata.axonId = repelem( ...
    data.summaryIds, cat(1, data.summary.nrChiasmata));
chiasmata.centerId = cat(1, data.summary.centerIdx);
chiasmata.summary = cat(1, data.summary.tracings);

chiasmata.solved = cat(1, data.summary.solved);
chiasmata.valid = cellfun(@(s) s.partitionIsValid, chiasmata.summary);

% restrict to unsolved chiasmata
chiasmata(chiasmata.solved, :) = [];
chiasmata(chiasmata.valid, :) = [];

%% export examples
rng(0);
rows = randperm(size(chiasmata, 1));

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for curIdx = 5
    curRow = rows(curIdx);
    curChi = chiasmata(curRow, :);
    curSummary = curChi.summary{1};
    
    % build original axon
    curAxonOld = data.oldAxons.axons(curChi.axonId);
    
    curComments = cell(size(curAxonOld.nodes, 1), 1);
    curComments(:) = {''};
    curComments{curChi.centerId} = 'Unsolved chiasma';
    
    curSkel = skeleton();
    curSkel = curSkel.addTree( ...
        'Original axon', curAxonOld.nodes(:, 1:3), ...
        curAxonOld.edges, [], [], curComments);
    
    % build split axons
    curAxonsNew = data.axons(data.parentIds == curChi.axonId);
    
    for curAxonIdx = 1:numel(curAxonsNew)
        curAxon = curAxonsNew(curAxonIdx);
        
        curSkel = curSkel.addTree( ...
            sprintf('Split axon %d', curAxonIdx), ...
            curAxon.nodes(:, 1:3), curAxon.edges);
    end
    
    % add flight paths
    for curExitIdx = 1:numel(curSummary.taskIds)
        curOverlaps = curSummary.overlaps{curExitIdx};
        curFlightName = sprintf('Flight (%d → %d)', curOverlaps);
        curFlightNodes = curSummary.nodes{curExitIdx};
        
        curSkel = curSkel.addTree(curFlightName, curFlightNodes);
    end
    
    curSkelName = sprintf( ...
        '%03d_axon-%d_node-%d.nml', ...
        curIdx, curChi.axonId, curChi.centerId);
    
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    % curSkel.write(fullfile(outputDir, curSkelName));
end