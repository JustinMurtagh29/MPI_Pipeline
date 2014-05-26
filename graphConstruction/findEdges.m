function [edges] = findEdges(nSegId)
% CAREFUL: assumes watershed based segmentation!!!
% e.g. according to 26 connectivity border width of 1 voxel set to 0, border voxel
% always has a neighbour that is wall voxel as well ....

% Alternate approach for finding unique neighbours for each wall voxel
nSegId = nSegId'; % should switch dimensionality in construction already if using this approach
nSegId = sort(nSegId,1);
lSegId = [false(1,size(nSegId,2)); diff(nSegId,1,1)>0];
nSegId = nSegId(lSegId);
list1 = mat2cell(nSegId, sum(lSegId,1), 1);

% Extract all pairs for higher level indices
fun2 = @(x) combnk(x,2);
pairs = cellfun(fun2,list1,'UniformOutput', false);
pairs = cell2mat(pairs);

% Remove duplicate pairs to find edges
edges = unique(pairs, 'rows');

end
