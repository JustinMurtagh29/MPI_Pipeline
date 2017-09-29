% filterPathLength 
% filters all agglomerates according to path length (designed for L4)
% -------------------------------------------------------------------------
% Inputs:
% agglos - agglomerates as sets of segments
% points - agglomerates as points from segment meta
% varThr - threshold for the variance explained by first PC
% spThr - threshold for the angle between the mean axis and the considered
% Outputs:
% filtered agglos & points, sps - scalarproducts to mean FPC
% meanFPC - mean first principal component
% rInd - indeces of filtered
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
function [ filteredAgglos, filteredPoints, meanFPC, sps, rInd ] = filterMainAxis( agglos, points, varThr, spThr, rInd )
tic;
% calculate PCA and var explained for each agglo
pcaMain = zeros(3, length(points)); 
varExplained = zeros(length(points),3);
for a=1:length(points)
    D = points{a};
    [coeff,score,latent,tsquared,explained,mu] = pca(D);
    pcaMain(:,a) = coeff(:,1);
    varExplained(a,:) = explained(1);
end
varExpFirst = varExplained(:,1);
disp('finished applying PCA to data.');

% kick out all agglos that are below threshold (don't explain enough
% variance) and get mean first principal component
varThr = varThr * 100;
varMask = varExpFirst >= varThr;
pcaMain = pcaMain(:,transpose(varMask));
filteredAgglos = agglos(varExpFirst >= varThr);
filteredPoints = points(varExpFirst >= varThr);
meanFPC = mean(pcaMain, 2);
rInd = rInd(varExpFirst >= varThr);

% filter that don't align with main axis
% use scalar product (angle)
sps = transpose(transpose(meanFPC) * pcaMain);
filteredAgglos = filteredAgglos(sps >= spThr);
filteredPoints = filteredPoints(sps >= spThr);
rInd = rInd(sps >= spThr);

f1 = length(filteredAgglos) / length(agglos);
fprintf('filtered agglos, fraction of total: %f (%d/%d)\n', f1, length(filteredAgglos), length(agglos));

disp('finished execution of filterMainAxis.');
toc; disp('---------------------------------------');


end

