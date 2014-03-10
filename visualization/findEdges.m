function edges = findEdges(seg)

ids = 1:max(seg(:));
edges = cell(length(ids),1);

% Structuring element for 6 connectivity between segments crossing borders
se = zeros(5,5,5);
se(2:4,2:4,2:4) = 1;
se(1,3,3) = 1;
se(5,3,3) = 1;
se(3,1,3) = 1;
se(3,5,3) = 1;
se(3,3,5) = 1;
se(3,3,1) = 1; 

for i=1:length(ids)
    % Find all neighbours
    object = seg == ids(i); 
    voxel = imdilate(object, se) & ~object;
    values = reshape(seg(voxel), sum(voxel(:)), 1, 1);
    uniqueValues = unique(values);
    uniqueValues(uniqueValues == 0) = [];
    edges{i} = zeros(length(uniqueValues), 2);
    for j=1:length(uniqueValues)
        edges{i}(j,:) = [ids(i) uniqueValues(j)];
    end
end
edges = cell2mat(edges);
edges = unique(edges, 'rows');
end


