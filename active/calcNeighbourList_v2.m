function g = calcNeighbourList_v2( g, raw, aff, seg)

ids = unique(seg);
ids(ids == 0) = [];
edges = cell(length(ids),1);
weights = cell(length(ids),1);
%matlabpool 12;
for i=1:length(ids)
    % Fing all neighbours
    object = seg == ids(i);
    voxel1 = bwdist(object, 'chessboard') < 3 & ~object;
    values = reshape(seg(voxel1), sum(voxel1(:)), 1, 1);
    uniqueValues = unique(values);
    uniqueValues(uniqueValues == 0) = [];
    % Find border of object
    voxel1 = bwdist(object, 'chessboard') < 2 & ~object;
    edges{i} = zeros(length(uniqueValues), 2);
    weights{i} = cell(length(uniqueValues), 1);
    % Loop over all neighbours and calculate some features for each edge
    for j=1:length(uniqueValues)
        edges{i}(j,:) = [ids(i) uniqueValues(j)];
        otherObject = seg == uniqueValues(j);
        voxel2 = bwdist(otherObject, 'chessboard') < 2 & ~otherObject;
        % Feature 1: Mean of affinity within border
        weights{i}{j,1} = aff(voxel1 & voxel2);
        % Feature 2: Surface between objects compared to 
        weights{i}{j,2} = sum(voxel1 & voxel2)/sum(object +otherObject);
    end
    visualizeObjectNeighbourhoodOnTheFly(raw, seg, edges{i}(:,2), weights{i}, ids(i));
end
%matlabpool close

for i=1:length(edges)
    g = g.insertEdges(edges{i}, weights{i});
end

end

