% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

resultsFile = '/home/amotta/Desktop/20170929T113827-results.mat';
outputDir = '/home/amotta/Desktop/chiasmata-splitting';
data = load(resultsFile);

%% Show fraction solved per size
clearvars -except outputDir data;

nrExits = cat(1, data.summary.nrExits);
solved = cat(1, data.summary.solved);

results = table;
[results.nrExits, ~, nrExits] = unique(nrExits);
results.nrSolved = accumarray(nrExits, solved, [], @sum, 0);
results.nrTotal = accumarray(nrExits, solved, [], @numel, 0);

disp(results)

%% Export examples
clearvars -except outputDir data;

nrExits = 4;
solved = false;

mkdir(outputDir);

prefix = {'unsolved', 'solved'};
prefix = sprintf('%d-fold__%s', nrExits, prefix{1 + solved});

t = table;
t.nrExits = cat(1, data.summary.nrExits);
t.solved = cat(1, data.summary.solved);
t.centerIdx = cat(1, data.summary.centerIdx);

allNrChiasmata = cat(1, data.summary.nrChiasmata);
t.axonId = repelem((1:numel(allNrChiasmata))', allNrChiasmata);
t.chiasmaId = cell2mat(arrayfun(@(ids) ...
    (1:ids)', allNrChiasmata, 'UniformOutput', false));

% find the ones which match criteria
t = t(t.nrExits == nrExits, :);
t = t(t.solved == solved, :);

% select random examples
rng(0);
t = t(randperm(size(t, 1), 10), :);

parentIds = find(data.oldAgglos.indBigAxons);
t.parentId = parentIds(t.axonId);

for row = 1:size(t, 1)
    skel = skeleton();
    
    % add parent
    parentId = t.parentId(row);
    agglo = data.oldAgglos.axons(parentId);
    
    comments = repmat({''}, size(agglo.nodes, 1), 1);
    comments{t.centerIdx(row)} = 'Chiasma';
    skel = skel.addTree( ...
        'Original', agglo.nodes, agglo.edges, [], [], comments);
    
    if solved
        % in "solved" case, show components
        newAgglos = data.newAgglos(data.parentIds == parentId);
        for newAgglo = newAgglos(:)'
            skel = skel.addTree('Split', newAgglo.nodes, newAgglo.edges);
        end
    else
        % in "unsolved" case, also show flight paths
        axonId = t.axonId(row);
        chiasmaId = t.chiasmaId(row);
        flights = data.summary(axonId).tracings{chiasmaId};
        
        for ffIdx = 1:numel(flights.nodes)
            ffNodes = flights.nodes{ffIdx};
            ffProcessed = flights.processed(ffIdx);
            
            ffOverlaps = flights.overlaps{ffIdx};
            ffOverlaps = arrayfun(@num2str, ffOverlaps, 'Uni', false);
            ffOverlaps = cat(2, ffOverlaps(:)', {'?', '?'});
            ffOverlaps = strjoin(ffOverlaps(1:2), '-');
            
            ffEdges = Graph.getMST(bsxfun( ...
                @times, ffNodes, data.p.raw.voxelSize));
            
            ffName = sprintf( ...
                'Flight Path %d. Processed %d. Overlaps %s', ...
                ffIdx, ffProcessed, ffOverlaps);
            skel = skel.addTree(ffName, ffNodes, ffEdges);
        end
    end
    
    skelFile = sprintf( ...
        '%s__axon-%d__chiasma-%d.nml', ...
        prefix, t.axonId(row), t.chiasmaId(row));
    skelFile = fullfile(outputDir, skelFile);
    
    skel = Skeleton.setParams4Pipeline(skel, data.p);
    skel.write(skelFile);
end
