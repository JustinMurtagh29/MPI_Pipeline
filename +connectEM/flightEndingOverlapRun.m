% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

%% load all the input data
m = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
param = m.p;

m = load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat', 'axonsNew');
origAgglos = m.axonsNew;

endings = load('/tmpscratch/scchr/AxonEndings/queriesBasedOnClusters/clusterData.mat');

m = load('/tmpscratch/scchr/AxonEndings/axonQueryResults/ff_struct_CS_MB_L4_AxonLeftQueries.mat', 'ff');
flights = m.ff;

m = load('/tmpscratch/scchr/AxonEndings/axonQueryResults/results_round2.mat');
flightResults = m.results;

m = load('/tmpscratch/kboerg/chiasmarunAugust/superagglos_postsplit.mat', 'superagglos');
superAgglos = m.superagglos;

clear m;

%% run main function
flightEndingOverlap = connectEM.flightEndingOverlap( ...
        param, origAgglos, endings, flights, flightResults, superAgglos);
    
%% debugging
%{
rng(0);
curRandId = find(flightEndingOverlap > 0);
curRandId = curRandId(randperm(numel(curRandId), 1));

curSuperAggloId = max(flightResults.endAgglo{curRandId});
curSuperAgglo = superAgglos(curSuperAggloId);

curEndingId = flightEndingOverlap(curRandId);
curEndingAggloIds = cat(1, 0, cumsum(cellfun(@max, endings.T)));
curEndingAggloIdx = find(curEndingAggloIds < curEndingId, 1, 'last');
curEndingIdx = curEndingId - curEndingAggloIds(curEndingAggloIdx);

curEndings = endings.borderPositions{curEndingAggloIdx};
curEndings = curEndings(endings.T{curEndingAggloIdx} == curEndingIdx, :);

curFlightNodes = flights.nodes{curRandId};

% build skeleton
curSkel = skeleton();
curSkel = curSkel.addTree( ...
    'Super-Agglomerate', curSuperAgglo.nodes(:, 1:3), curSuperAgglo.edges);
curSkel = curSkel.addTree( ...
    sprintf('Flight #%d', curRandId), curFlightNodes);

% add endings
curSkel = curSkel.addNodesAsTrees(curEndings);
curSkel.names((end - size(curEndings, 1) + 1):end) = {'Ending'};

curFileName = sprintf('%d_flight.nml', curRandId);
curFileName = fullfile('/home/amotta/Desktop', curFileName);

curSkel = Skeleton.setParams4Pipeline(curSkel, param);
curSkel.write(curFileName);
%}