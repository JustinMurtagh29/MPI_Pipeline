% extractApicals 
% script to filter a state of dendrite agglos for apicals (designed for L4)
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>


% load necessary data for filtering
dendAggloStatePath = '/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/dendrites_03.mat';
spineHeadsStatePath = '/gaba/tmpscratch/zecevicm/L4/Apicals/20170829_spineheadsAttached.mat';
dendLensPath = '/gaba/tmpscratch/zecevicm/L4/Apicals/20170829_DendriteAggloLengths.mat';
parameterPath = '/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat';
metaPath = '/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat';
local = '';
[ agglosNew, agglos, agglos_all, points, indBigDends, bounds, spineHeads, dendLens, sM, p ] = ...
    Apicals.loadDendAggloState(dendAggloStatePath, spineHeadsStatePath, dendLensPath, parameterPath, metaPath, local);


% filter using site criterion
tol = 200;
[ filteredAgglos, filteredPoints, rInd_1 ] = Apicals.filterSiteCriterion( agglos, points, bounds, tol );


% filter using path length
ulbound = [200000 50000];
[ filteredAgglos, filteredPoints, pathLengths, rInd_2 ] = ...
                                            Apicals.filterPathLength( filteredAgglos, filteredPoints, ulbound, rInd_1 );
% plot: histogram for pathLength distribution
h = histogram(pathLengths);
title([ num2str(length(pathLengths)) ' agglos out of ' num2str(length(agglos)) ' - Bin width: ' num2str(h.BinWidth)]);
xlabel('physical path length'); ylabel('count');


% (two steps) filter using PCA (variance explained) and main axis alignment
varThr = 0.98;
spThr = 0.95;
[ filteredAgglos, filteredPoints, sps, rInd_3 ] = ...
                                        Apicals.filterMainAxis( filteredAgglos, filteredPoints, varThr, spThr, rInd_2 );
% plot: histogram for distribution of distances (meanFPC from observed)
h = histogram(sps);
title([ num2str(length(sps)) ' agglos out of ' num2str(length(agglos)) ' - Bin width: ' num2str(h.BinWidth)]);
xlabel('distance from mean FPC'); ylabel('count');


% collect spine head counts for all agglos (!) because of attachment
% seems like the agglos are sorted the way that > 5um come first
% would explain why below worked with this as before
[ shCounts ] = Apicals.collectSHcounts( agglos_all, spineHeads );

% collect spine head densities for the filtered Agglos
shDensity = shCounts(rInd_3) ./ dendLens(rInd_3);
% replace NaN with zero spine density
shDensity(isnan(shDensity)) = 0; % zero spines
% plot: histogram for spine head densities
histogram(shDensity, 10);
xlabel('Spine density: # spines per um'); ylabel('# dendrite agglos');
title(['Spine Head density for ' num2str(length(rInd_3)) ' agglos filtered']);

