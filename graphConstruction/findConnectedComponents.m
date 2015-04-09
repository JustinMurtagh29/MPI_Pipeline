function [equivalenceClasses, objectClassLabels] = findConnectedComponents(sparseAdjMatrixOrEdgeList)
% Finds connected components in undirected graph (graph is made symetric)
%INPUT: either sparse adjanceny matrix (assumed symetric) or one direction edge list
% (which will be made into symetric adjanceny list here, inconsistent?)
%OUTPUT: first result = cell array with all equivalence classes with more
% than 1 object, second result equivalence class label for each object

% Create sparse undirected adjaceny matrix if not supplied
if ~issparse(sparseAdjMatrixOrEdgeList)
    sparseAdjMatrixOrEdgeList = double(sparseAdjMatrixOrEdgeList); % integer sparse matrices apparently not possible in matlab (tested R2011b)
    maxValue = max(sparseAdjMatrixOrEdgeList(:));
    sparseAdjMatrixOrEdgeList = sparse(sparseAdjMatrixOrEdgeList(:,1),sparseAdjMatrixOrEdgeList(:,2),1,maxValue,maxValue);
end

% Make sure adjMatric is symetric
sparseAdjMatrixOrEdgeList = sparseAdjMatrixOrEdgeList + sparseAdjMatrixOrEdgeList';

% Find block diagonal matrix permutation
[rowPermutation,~,rowBlockBoundaries] = dmperm(sparseAdjMatrixOrEdgeList + speye(maxValue));

% Create vector with one at each row boundary start
newLabelStart = zeros(1,maxValue);
newLabelStart(rowBlockBoundaries(1:end-1)) = 1;

% Calculate object class labels (equivalence class of each object)
objectClassLabels = cumsum(newLabelStart);
objectClassLabels(rowPermutation) = objectClassLabels;

% Equivalence classes (sort out 1-element ones first)
sizeBlocks = diff(rowBlockBoundaries);
sizeBlocks(sizeBlocks == 1) = [];
rowPermutation(~(newLabelStart == 0 | [diff(newLabelStart) 0] == -1)) = [];
equivalenceClasses = mat2cell(rowPermutation', sizeBlocks);

end
