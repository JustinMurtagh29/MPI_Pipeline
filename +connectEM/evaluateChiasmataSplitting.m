% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

resultsFile = '/home/amotta/Desktop/20170928T222837-results.mat';
data = load(resultsFile);

%% Show fraction solved per size
clearvars -except data;

nrExits = cat(1, data.summary.nrExits);
solved = cat(1, data.summary.solved);

results = table;
[results.nrExits, ~, nrExits] = unique(nrExits);
results.nrSolved = accumarray(nrExits, solved, [], @sum, 0);
results.nrTotal = accumarray(nrExits, solved, [], @numel, 0);

disp(results)

%% Export examples
clearvars -except data;

nrExits = 5;
solved = true;
outputDir = '/home/amotta/Desktop';

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
    
    % add components
    newAgglos = data.newAgglos(data.parentIds == parentId);
    for newAgglo = newAgglos(:)'
        skel = skel.addTree('Split', newAgglo.nodes, newAgglo.edges);
    end
    
    skelFile = sprintf( ...
        '%s__axon-%d__chiasma-%d.nml', ...
        prefix, t.axonId(row), t.chiasmaId(row));
    skelFile = fullfile(outputDir, skelFile);
    
    skel = Skeleton.setParams4Pipeline(skel, data.p);
    skel.write(skelFile);
end