function g = calcNeighbourList( g, raw, aff, seg)

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
    % Loop over all neighbours and save histogram of affinity on overlap
    for j=1:length(uniqueValues)
        edges{i}(j,:) = [ids(i) uniqueValues(j)];
        otherObject = seg == uniqueValues(j);
        voxel2 = bwdist(otherObject, 'chessboard') < 2 & ~otherObject;
        weights{i}{j} = aff(voxel1 & voxel2);
    end
    visualizeObjectNeighbourhoodOnTheFly(raw, seg, edges{i}(:,2), weights{i}, ids(i));
end
%matlabpool close

for i=1:length(edges)
    g = g.insertEdges(edges{i}, weights{i});
end

end

