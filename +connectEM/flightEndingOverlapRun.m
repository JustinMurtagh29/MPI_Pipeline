function flightEndingOverlap = flightEndingOverlapRun()

% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

%% load all the input data
m = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
param = m.p;

% m = load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat', 'axonsNew');
% origAgglos = m.axonsNew;
m = load(fullfile(param.saveFolder, 'aggloState/', 'axons_04.mat'));
origAgglos = m.axons;

endings = load(fullfile(param.saveFolder, 'aggloState/', 'axonEndings.mat'));

% m = load('/tmpscratch/scchr/AxonEndings/axonQueryResults/ff_struct_CS_MB_L4_AxonLeftQueries.mat', 'ff');
m = load(fullfile(param.saveFolder, 'aggloState/', 'AxonFlightPaths.mat'), 'ff');
flights = m.ff;

% m = load('/tmpscratch/scchr/AxonEndings/axonQueryResults/results_round2.mat');
m = load(fullfile(param.saveFolder, 'aggloState/', 'AxonQueryOverlaps.mat'), 'results')
flightResults = m.results;

% m = load('/tmpscratch/kboerg/chiasmarunAugust/superagglos_postsplit.mat', 'superagglos');
m = load(fullfile(param.saveFolder, 'aggloState/', 'axons_04.mat'));
superAgglos = m.axons;

clear m;

%% run main function
flightNodes = flights.nodes;
flightAgglos = cellfun( ...
    @union, flightResults.startAgglo, ...
    flightResults.endAgglo, 'UniformOutput', false);

flightEndingOverlap = connectEM.flightEndingOverlap( ...
        param, origAgglos, endings, flightNodes, flightAgglos, superAgglos);
