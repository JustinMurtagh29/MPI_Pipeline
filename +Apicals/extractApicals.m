% extractApicals 
% extracts apical dendrite candidates from agglomerates (designed for L4)
% -------------------------------------------------------------------------
% Inputs:
% dendAggloStatePath - current state of dendrite agglomerates
% spineHeadsStatePath - current state of attached spine heads
% dendLensPath - lengths for the agglomerates (in um)
% parameterPath - parameter file
% metaPath - segment meta file
% local - mount point when working locally e.g. '/mnt/whatever'
% tol - tolerance on the contact site filtering (in pixel, e.g. 100~1um)
% ulbound - 2x1 upperbound for length and lower bound (in pixel)
% varThr - threshold for variance explained
% spThr - threshold for scalarproduct (main axis alignment)
% Outputs:
% filtered agglos, remainingIds - indices wrt file, spineHeads, dendLens -
% lengths of dendrites (in um)
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
function [filteredAgglos, remainingIds, spineHeads, dendLens] = ...
    extractApicals(dendAggloStatePath, spineHeadsStatePath, dendLensPath, parameterPath, metaPath, local,...
                    tol, ulbound, varThr, spThr)

% load necessary data for filtering
[ agglosNew, agglos, agglos_all, points, indBigDends, bounds, spineHeads, dendLens, sM, p ] = ...
    Apicals.loadDendAggloState(dendAggloStatePath, spineHeadsStatePath, dendLensPath, parameterPath, metaPath, local);


% filter using site criterion
[ filteredAgglos, filteredPoints, rInd_1 ] = Apicals.filterSiteCriterion( agglos, points, bounds, tol );


% filter using path length
[ filteredAgglos, filteredPoints, pathLengths, rInd_2 ] = ...
                                            Apicals.filterPathLength( filteredAgglos, filteredPoints, ulbound, rInd_1 );


% (two steps) filter using PCA (variance explained) and main axis alignment
[ filteredAgglos, filteredPoints, sps, rInd_3 ] = ...
                                        Apicals.filterMainAxis( filteredAgglos, filteredPoints, varThr, spThr, rInd_2 );

% collect indeces of the agglos after each filtering step
remainingIds = {rInd_1, rInd_2, rInd_3};
                                    
% filtering with spine head density possibly to be implemented later on

end
