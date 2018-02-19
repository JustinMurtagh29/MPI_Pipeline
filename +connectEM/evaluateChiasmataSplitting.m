% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
resultsFile = fullfile(rootDir, 'chiasmataSplitting/20171009T193744-kmb-on-axons-6c/outputs/20171030T181930_results.mat');
outputDir = '/home/amotta/Desktop/chiasmata-splitting';

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

data = load(resultsFile);

mkdir(outputDir);

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
solved = true;

prefix = {'unsolved', 'solved'};
prefix = sprintf('%d-fold__%s', nrExits, prefix{1 + solved});

t = table;
t.nrExits = cat(1, data.summary.nrExits);
t.solved = cat(1, data.summary.solved);
t.centerIdx = cat(1, data.summary.centerIdx);

allNrChiasmata = cat(1, data.summary.nrChiasmata);
t.axonId = repelem(data.summaryIds, allNrChiasmata);
t.chiasmaId = cell2mat(arrayfun(@(ids) ...
    (1:ids)', allNrChiasmata, 'UniformOutput', false));

% find the ones which match criteria
t = t(t.nrExits == nrExits, :);
t = t(t.solved == solved, :);

% select random examples
rng(0);
t = t(randperm(size(t, 1), 100), :);

parentIds = find(data.oldAxons.indBigAxons);
t.parentId = parentIds(t.axonId);

for row = 1:size(t, 1)
    skel = skeleton();
    
    % add parent
    parentId = t.parentId(row);
    agglo = data.oldAxons.axons(parentId);
    
    comments = repmat({''}, size(agglo.nodes, 1), 1);
    comments{t.centerIdx(row)} = 'Chiasma';
    skel = skel.addTree( ...
        'Original', agglo.nodes, agglo.edges, [], [], comments);
    
    if solved
        % in "solved" case, show components
        newAgglos = data.axons(data.parentIds == parentId);
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
        '%0*d_%s__axon-%d__chiasma-%d.nml', ...
        ceil(log10(1 + size(t, 1))), row, ...
        prefix, t.axonId(row), t.chiasmaId(row));
    skelFile = fullfile(outputDir, skelFile);
    
    skel = Skeleton.setParams4Pipeline(skel, data.p);
    skel.write(skelFile);
end

%% Compare size of largest agglomerates (in nodes)
clearvars -except outputDir param data;

buildAggloSizes = @(axons) sort(arrayfun( ...
    @(a) size(a.nodes, 1), axons), 'descend');
aggloSizesNew = buildAggloSizes(data.axons);
aggloSizesOld = buildAggloSizes(data.oldAxons.axons);

largestAggloSizes = table;
largestAggloSizes.before = aggloSizesOld(1:20);
largestAggloSizes.after = aggloSizesNew(1:20);

fprintf('\n');
fprintf('Largest agglomerates (in nodes) before / after chiasma splitting:\n\n');
disp(largestAggloSizes);

%% Export largest remaining agglomerates (in nodes)
clearvars -except outputDir param data;

agglos = data.axons;
aggloSizes = arrayfun(@(a) size(a.nodes, 1), agglos);

[~, sortIds] = sort(aggloSizes, 'descend');
agglos = agglos(sortIds);

for curIdx = 1:10
    curAggloId = sortIds(curIdx);
    curAgglo = agglos(curIdx);
    
    skel = skeleton();
    skel = skel.addTree( ...
        sprintf('Agglomerate #%d', curAggloId), ...
        curAgglo.nodes, curAgglo.edges);
    skel = Skeleton.setParams4Pipeline(skel, param);
    
    skelFile = sprintf('largest-%d__axon-%d.nml', curIdx, curAggloId);
    skel.write(fullfile(outputDir, skelFile));
end
