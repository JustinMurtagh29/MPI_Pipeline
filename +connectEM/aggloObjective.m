function value = aggloObjective( probMatrix, edgeMatrix )
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
intraProbMatrix = ccSortedProbMatrix(ccBlockMatrix);

intraProbability = sum(nonzeros(intraProbMatrix));
intraProbability
interProbability = sum(nonzeros(ccSortedProbMatrix)) - intraProbability;
interProbability
value = interProbability - intraProbability;

% Get nonzero elements and positions for summing (tried indexing in sparse
% matrices but seems to be very costly, at least more costly than this)
% [r, c, v] = find(ccSortedProbMatrix);
% withinComponent = false(size(r));
% for i=1:length(r)
%     idx1 = r(i) > blockStart & r(i) < blockEnd;
%     if any(idx1)
%         idx2 = c(i) > blockStart & c(i) < blockEnd;
%         if any(idx2) && all(idx1 == idx2)
%             withinComponent(i) = true;
%         end
%     end
% end

% % Calculate intra und inter component probability sum
% intraProbability = zeros(length(rowBlockBoundaries)-1,1);
% for i=1:length(rowBlockBoundaries)-1
%     if sizeBlocks(i) > 1
%         display(num2str(sizeBlocks(i)));
%         tic;
%         idxStart = rowBlockBoundaries(i);
%         idxEnd = rowBlockBoundaries(i+1)-1;
%         tempMatrix = ccSortedProbMatrix(idxStart:idxEnd, idxStart:idxEnd);
%         intraProbability(i) = sum(tempMatrix(:));
%         toc;
%     end
% end
% 
% % Something like this is what we might want as a metric?
% intraProbability = sum(intraProbability);
% interProbability = sum(ccSortedProbMatrix(:)) - intraProbability;
% value = interProbability - intraProbability;

end
