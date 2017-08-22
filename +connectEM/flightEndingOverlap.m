% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

%% load dataset parameters and super-agglomerates
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat', 'p');
load('/tmpscratch/kboerg/chiasmarunAugust/superagglos_postsplit.mat', 'superagglos');

%% build agglomerates
agglos = arrayfun(@(sa) sa.nodes(:, 4), superagglos, 'Uni', false);
agglos = cellfun(@(ids) unique(ids(~isnan(ids))), agglos, 'Uni', false);

%% load flight paths and their evaluation
load('/tmpscratch/scchr/AxonEndings/axonQueryResults/ff_struct_CS_MB_L4_AxonLeftQueries.mat', 'ff');
data = load('/tmpscratch/scchr/AxonEndings/axonQueryResults/results_round2.mat');

%% group flight paths per agglomerate
endAgglos = data.results.endAgglo;
endAgglos = cellfun( ...
    @(a) reshape(a, [], 1), endAgglos, 'Uni', false);
flightIds = repelem( ...
    reshape(1:numel(endAgglos), [], 1), cellfun(@numel, endAgglos));

aggloFlightIds = accumarray( ...
    cell2mat(endAgglos), flightIds, [], @(ids) {ids(:)}, {});
aggloWithFlighIds = find(~cellfun(@isempty, aggloFlightIds));

%% load endings
% these are the original axon agglomerates
load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat', 'axonsNew');

% build look-up table
maxSegId = Seg.Global.getMaxSegId(p);
axonsNewLUT = zeros(maxSegId, 1);
axonsNewLUT(cell2mat(axonsNew)) = repelem( ...
    1:numel(axonsNew), cellfun(@numel, axonsNew));

% load border CoMs
load(fullfile(p.saveFolder, 'globalBorder.mat'), 'borderCoM');

% load endings
data = load('/tmpscratch/scchr/AxonEndings/queriesBasedOnClusters/clusterData.mat');
axonNewBorderIds = cell(numel(axonsNew), 1);
axonNewBorderIds(data.mapping) = data.thisBorderIdx;

%% generate
rng(13);
curAggloId = randperm(numel(aggloWithFlighIds), 1);
curAggloId = aggloWithFlighIds(curAggloId);

curSuperAgglo = superagglos(curAggloId);
curSegIds = curSuperAgglo.nodes(:, 4);
curSegIds(isnan(curSegIds)) = [];
curSegIds = setdiff(curSegIds, 0);

curOrigAgglos = setdiff(axonsNewLUT(curSegIds), 0);
curEndingBorderIds = axonNewBorderIds(curOrigAgglos(:));
curEndingBorderIds = cell2mat(curEndingBorderIds);

curSkel = skeleton();
curSkel = curSkel.addTree( ...
    'Super-Agglomerate', curSuperAgglo.nodes(:, 1:3), curSuperAgglo.edges);

% add flights
for curFlightId = reshape(aggloFlightIds{curAggloId}, 1, [])
    curName = sprintf('Flight #%d', curFlightId);
    curSkel = curSkel.addTree(curName, ff.nodes{curFlightId});
end

% add endings
for curBorderIdx = reshape(curEndingBorderIds, 1, [])
    curName = sprintf('Ending (Border #%d)', curBorderIdx);
    curPos = borderCoM(curBorderIdx, :);
    curSkel = curSkel.addTree(curName, curPos);
end

curFileName = sprintf('%d_superagglo-with-inbound.nml', curAggloId);
curFileName = fullfile('/home/amotta/Desktop', curFileName);

curSkel = Skeleton.setParams4Pipeline(curSkel, p);
curSkel.write(curFileName);