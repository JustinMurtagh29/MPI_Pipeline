function calculateSurfaceInInnerCube( skelPath, skelFile, outputFile )
    display(['Processing skeleton: ' skelFile]);
    [~,skel_data] = evalc('parseNml([skelPath skelFile])');
    nodes = skel_data{1,1}.nodes(:,1:3);
    % for each node, find cube in which it lays so the cubes data can be used for
    % several nodes
    nodeData.nodes = nodes;
    nodeData.cubeCoords = zeros(1,size(nodes,2));
    for i = 1 : size(nodes,1)
        nodeData.cubeCoords(i,:) = floor(( nodes(i,:) - 1) / 128 ); %from readKnossosRoi ...  overlap of 128!
        nodeData.flag(i) = 0;
    end
    %group nodes that lie in the same cube
    groupCount = 0;
    for i = 1 : size(nodeData.nodes,1)
        if(nodeData.flag(i)==0)
            groupCount = groupCount+1;
            groupedNodes{1,groupCount}.nodes = nodeData.nodes(i,:);
            groupedNodes{1,groupCount}.cubeCoords = nodeData.cubeCoords(i,:);
            if(i<size(nodeData.nodes,1))            
                for j = i+1 : size(nodeData.nodes,1)
                    if(nodeData.cubeCoords(i,:) == nodeData.cubeCoords(j,:))
                        groupedNodes{1,groupCount}.nodes = [groupedNodes{1,groupCount}.nodes ; nodeData.nodes(j,:)];
                        nodeData.flag(j) = 1;
                    end
                end
            end
        end
    end

    % for each cube get data and write 1 into result cube for the voxels 
    % belonging to one of the segments of the nodes lying in the cube
    area = zeros(29, 38, 41);
    for i = 1 : size(groupedNodes,2)
        %read cube
        if all(groupedNodes{i}.cubeCoords > [7 3 1]) & all(groupedNodes{i}.cubeCoords < [30 39 42])
            cube = readKnossosCube('/nfs/bmo/mberning/20140310backup/mag1/', '100527_k0563_seg', groupedNodes{i}.cubeCoords, 'uint16', '', 'raw', 256);
            %get the color values of the nodes in the cube
            segIds = zeros(1,size(groupedNodes{i}.nodes,1));
            zeroOfCube = groupedNodes{i}.cubeCoords * 128 + 1;
            zeroOfCube = zeroOfCube-64;
            for j = 1 : size(groupedNodes{i}.nodes,1) 
                rel_coords = groupedNodes{i}.nodes(j,:) - zeroOfCube(1,:); % point in cube = actual point - zero of cube
                segIds(1,j) = cube(rel_coords(1),rel_coords(2),rel_coords(3)); 
            end   
            segIds(segIds == 0) = []; %delete the zeros no neuron has color black
            counts = histc(segIds,unique(segIds));
            segIds = unique(segIds);
            segIds(counts < 1) = [];
            if ~isempty(segIds)
                for k = 1 : size(segIds,2)
                    cube(cube == segIds(k)) = NaN;
                end
                cube(~isnan(cube)) = 0;
                cube(isnan(cube)) = 1;
                cube = imclose(cube, ones([3,3,3]));
                cube = padarray(cube, [1 1 1]);
                cube = smooth3(cube, 'gaussian', 5, 2);
                cube = cube(65:end-64,65:end-64,65:end-64);
                issfs{i} = isosurface(cube, .5);
                if ~isempty(issfs{i}.vertices)
                    vertices = bsxfun(@times, issfs{i}.vertices, [12 12 25]);
                    faces = issfs{i}.faces;
                    a = vertices(faces(:, 2), :) - vertices(faces(:, 1), :);
                    b = vertices(faces(:, 3), :) - vertices(faces(:, 1), :);
                    c = cross(a, b, 2);
                    cubeCoord = groupedNodes{i}.cubeCoords;
                    area(cubeCoord(1),cubeCoord(2),cubeCoord(3))= 1/2 * sum(sqrt(sum(c.^2, 2)));               
                end
            end
        end
    end
    save(outputFile, 'area');
end

