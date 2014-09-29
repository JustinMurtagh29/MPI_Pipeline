function output = seededReconstruction(adjMatrix, startSv)
%not finished
numNodesFirst = 5;
numNodesSecond = 5;

[B, idx] = sort(adjMatrix(startSv, :));
output = idx(1:numNodesFirst);

%find neighbors
for i = 1:numNodesFirst
    [B, idx] = sort(adjMatrix(output(i),:));
    neighbors(:,i) = idx(1:numNodesSecond); 
end

%calculate cumulative prob

end

