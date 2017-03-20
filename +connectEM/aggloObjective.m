function [value, edgeMatrixNew, intraCostPos, intraCostNeg, interCostPos, interCostNeg] = aggloObjective( probMatrix, edgeMatrix )
% Calculate probability within and between connected components
% Similar to Graph.findConnectedComponents but optimized for this specific
% task (and its speed)

% Find permutation needed for block diagonal form & block boundaries
[rowPermutation, ~, rowBlockBoundaries, ~] = dmperm(edgeMatrix + ...
    edgeMatrix' + speye(size(edgeMatrix)));
% Modify probMatrix into this block diagonal form, now probMatrix has
% rectangles along diagonal with components
ccSortedProbMatrix = probMatrix(rowPermutation, rowPermutation);

% Remove single segments from rowBlockBoundaries as they cannot have a
% probability with themselves
sizeBlocks = diff(rowBlockBoundaries);
idx = sizeBlocks > 1;
blockStart = rowBlockBoundaries(idx);
blockEnd = blockStart + sizeBlocks(idx) - 1;

% Construct logical sparse matrix with blocks
blockIdx = arrayfun(@(x,y)combnk(x:y,2), blockStart, blockEnd, 'uni', false);
blockIdx = cat(1, blockIdx{:});
ccBlockMatrix = sparse(cat(1, blockIdx(:,1), blockIdx(:,2)), ...
    cat(1, blockIdx(:,2), blockIdx(:,1)), true, ...
    size(ccSortedProbMatrix,1), size(ccSortedProbMatrix,2));
intraProbMatrix = ccSortedProbMatrix(and(ccSortedProbMatrix,ccBlockMatrix));

% Construct new edge matrix (given CC)
rowPermutationInv(rowPermutation) = 1:length(rowPermutation);
edgeMatrixNew = edgeMatrix;
% TODO
% ccBlockMatrixOld = and(ccBlockMatrix(rowPermutationInv, rowPermutationInv), probMatrix);
% edgeMatrixNew(ccBlockMatrixOld) = true;
% clear ccBlockMatrixOld;

% within components reward high probability edges and punish low
% probability edges, vice-versa between components
intraCosts = nonzeros(intraProbMatrix) - 0.5;
intraCostPos = sum(intraCosts(intraCosts > 0));
intraCostNeg = 100 * sum(intraCosts(intraCosts < 0));
interCosts = nonzeros(ccSortedProbMatrix) - 0.5;
interCostPos = sum(interCosts(interCosts > 0)) - intraCostPos;
interCostNeg = sum(interCosts(interCosts < 0)) - intraCostNeg;
% Idea behind chosing the 2 weights above: intraCostNeg should be really
% expensive (= below 50% probability within component), 1e7 to make all
% terms similar order of magnitude, not sure whether important
value = interCostPos + interCostNeg ...
    - intraCostPos - intraCostNeg;

end
