% cutOutSoma 
% cuts sphere around center of agglomerate to delete soma (designed for L4)
% -------------------------------------------------------------------------
% Inputs:
% sM - segment meta file
% soma - soma agglomerate, agglo - agglomerate cotaining the soma
% thr - threshold for cutting the sphere
% Outputs:
% agglo with sphere cut around soma, mask used to cut out soma
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
function [ aggloWithoutSoma, mask ] = cutOutSoma( sM, soma, agglo, thr )
tic;
% get nodes for agglos
somaP = sM.point(:,soma)';
aggloP = sM.point(:,agglo)';

% get center by mean of the node values
center = round([mean(somaP(:,1)) mean(somaP(:,2)) mean(somaP(:,3))]);

% get distances
dist = pdist2(center, aggloP);

% cutout  (should roughly be the difference of the lengths I guess)
mask = dist >= thr;
aggloWithoutSoma = agglo(mask);


disp('finished execution of cutOutSoma.');
toc; disp('---------------------------------------');

end

