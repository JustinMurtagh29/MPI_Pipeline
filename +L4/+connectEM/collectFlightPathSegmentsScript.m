%% add segments along flight path to axons
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%% load data created in connectEM.axonQueryAnalysis

outputFolder = '/tmpscratch/mberning/axonQueryResults/';
m = load([outputFolder 'postQueryState.mat']);

%% get new agglos

axonsNew = L4.connectEM.addFlightPathSegmentsToAgglos(m.axons, ...
    m.uniqueSegments, m.eqClassCCfull, m.startAgglo, m.idxGood);

%% write 10 random examples to nml

diffIdx = find(cellfun(@length, axonsNew) ~= cellfun(@length, m.axonsNew));
fprintf('%.2f of axons have changed.\n', length(diffIdx)/length(axonsNew));

% rng('shuffle');
% idx = randi(length(diffIdx), 10, 1);
% result from random index
idx = [1170;1256;3532;4866;760;6159;2134;3092;6072;7148];

m2 = load(['/gaba/u/mberning/results/pipeline/20170217_ROI/' ...
    'segmentMeta.mat'], 'point');
point = m2.point';
clear m2

exAgglos = [axonsNew(diffIdx(idx)), m.axonsNew(diffIdx(idx))]';
exAgglos = exAgglos(:);
skel = L4.Agglo.agglo2Nml(exAgglos, point');

for i = 1:skel.numTrees()/2
    skel.names{2*i -1} = sprintf('Axon%02d_queryOnly', i);
    skel.names{2*i} = sprintf('Axon%02d_flightPathSegAdded', i);
end
skel.write('FlightPathSegmentCollection.nml');
