% A script for comparing the two chiasmata detections written by Kevin M.
% Boergens and Manuel Berning, respectively. The former primarily relies on
% graph topology while the latter one exploits the geometry properties of
% nodes and edges.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

%% load input data
% data was loaded with +connectEM/evaluateChiasmataDetection.m

dataA = load('/home/amotta/Desktop/input-data-kmb.mat');
dataA = dataA.chiasmata;

dataB = load('/home/amotta/Desktop/input-data-mb.mat');
dataB = dataB.chiasmata;

assert(isequal(size(dataA), size(dataB)));

%% 
aggloCount = numel(dataA);

nodesVenn = table;
nodesVenn.onlyA = zeros(aggloCount, 1);
nodesVenn.aAndB = zeros(aggloCount, 1);
nodesVenn.onlyB = zeros(aggloCount, 1);

chiasmaVenn = table;
chiasmaVenn.onlyA = cell(aggloCount, 1);
chiasmaVenn.aAlsoB = cell(aggloCount, 1);
chiasmaVenn.onlyB = cell(aggloCount, 1);
chiasmaVenn.bAlsoA = cell(aggloCount, 1);

tic;
for aggloIdx = 1:aggloCount
    nodeCount = size(dataA{aggloIdx}.nodes, 1);
    luts = zeros(nodeCount, 2);
    
    maskA = dataA{aggloIdx}.isIntersection;
    maskB = dataB{aggloIdx}.isIntersection;
    
    nodesVenn.onlyA(aggloIdx) = sum( maskA & ~maskB);
    nodesVenn.aAndB(aggloIdx) = sum( maskA &  maskB);
    nodesVenn.onlyB(aggloIdx) = sum(~maskA &  maskB);
    
    nodesA = dataA{aggloIdx}.ccNodeIdx;
    nodesB = dataB{aggloIdx}.ccNodeIdx;
    
    if ~isempty(nodesA)
        luts(cell2mat(nodesA(:)), 1) = repelem( ...
            (1:numel(nodesA))', cellfun(@numel, nodesA));
    end
    
    if ~isempty(nodesB)
        luts(cell2mat(nodesB(:)), 2) = repelem( ...
            (1:numel(nodesB))', cellfun(@numel, nodesB));
    end
    
    luts = luts(all(luts, 2), :);
   [pairs, ~, pairCount] = unique(luts, 'rows');
    pairCount = accumarray(pairCount, 1);
    
    chiasmaVenn.aAlsoB{aggloIdx} = unique(pairs(:, 1));
    chiasmaVenn.onlyA{aggloIdx} = setdiff( ...
        (1:numel(nodesA))', chiasmaVenn.aAlsoB{aggloIdx});
    chiasmaVenn.bAlsoA{aggloIdx} = unique(pairs(:, 2));
    chiasmaVenn.onlyB{aggloIdx} = setdiff( ...
        (1:numel(nodesB))', chiasmaVenn.bAlsoA{aggloIdx});
end
toc;

%%
results = table;
results.nrOnlyA = sum(nodesVenn.onlyA);
results.nrAAndB = sum(nodesVenn.aAndB);
results.nrOnlyB = sum(nodesVenn.onlyB);

fprintf('\nResults in terms of chiasmatic nodes:\n\n');
disp(results);

results = table;
results.nrOnlyA = sum(cellfun(@numel, chiasmaVenn.onlyA));
results.nrAAlsoB = sum(cellfun(@numel, chiasmaVenn.aAlsoB));
results.nrOnlyB = sum(cellfun(@numel, chiasmaVenn.onlyB));
results.nrBAlsoA = sum(cellfun(@numel, chiasmaVenn.bAlsoA));

results.nrB = results.nrOnlyB + results.nrBAlsoA;
results.nrA = results.nrOnlyA + results.nrAAlsoB;
results = results(:, [end, 1:(end - 1)]);

fprintf('\nResults in terms of chiasmata:\n\n');
disp(results)