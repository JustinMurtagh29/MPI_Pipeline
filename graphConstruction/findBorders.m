function [newedges,borders] = findBorders(ind, nInd, nSegId, edges, M, N, P)

% Calculate border voxel indices to nSegId for each unique segment in edges
allSegments = unique(edges(:));
idx = cell(length(allSegments),1);
for i=1:length(allSegments)
    idx{i} = find(any(nSegId==allSegments(i),2));
end

% Preallocation
borders = struct('PixelIdxList', {}, 'Area', {}, 'Centroid', {});
newedges = zeros(size(edges,1), 2, 'uint32'); % Preallocate biggest chunk, will grow due to multiple borders between objects, better approach?
g = 1;

for i=1:size(edges,1)
    edgeIdx1 = edges(i,1) == allSegments;
    edgeIdx2 = edges(i,2) == allSegments;
    
    % Find all borderpixel from edges (i,1) to edges(i,2) and save them in borderIdx
    temp = intersect(idx{edgeIdx1},idx{edgeIdx2});
    borderIdx = ind(temp);
    
    % Adjaceny list construction & BFS
    adjacency_list = nInd(temp,:);
    adjacency_list = mat2cell(adjacency_list, ones(size(adjacency_list,1),1), 26);
    adjacency_list = cellfun(@(x)intersect(x,borderIdx), adjacency_list, 'UniformOutput', false);
    bordersForThisEdge = bfs(borderIdx,adjacency_list);

    for j=1:length(bordersForThisEdge)
        % Reshaping the indices of regionmatrix to indices of seg
        [x y z] = ind2sub([M N P],bordersForThisEdge{j});
        x = int32(x)-1; % -1 due to padding in parent function (bad practice?)
        y = int32(y)-1;
        z = int32(z)-1;
        % Collect output(s)
        borders(g).PixelIdxList = sub2ind([M-2,N-2,P-2],x,y,z);
        borders(g).Area = size(borders(g).PixelIdxList,2);
        borders(g).Centroid = mean(vertcat(x,y,z),2)';
        newedges(g,:) = edges(i,:);
        g=g+1;
    end

end

end
