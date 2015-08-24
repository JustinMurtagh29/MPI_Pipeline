function galleryCortex(p, skelPath, skelFile, outputPath)

    % Display status
    display(['Processing skeleton: ' skelFile]);
    % If 2nd argument is string use normal mode of getting passed nml, otherwise expect skeleton struct as 2nd input
    if isa(skelPath, 'char')
        % evalc to supress output of parseNml
        [~,skel_data] = evalc('parseNml([skelPath skelFile])');
    else
        skel_data = skelPath;
    end

    % Find X,Y,Z for each node if it is within parameter.local(X,Y,Z).bboxSmall
    for sk = 1:length(skel_data)
        if size(skel_data{sk}.nodes,1) > 10 
            nodeData.nodes = skel_data{sk}.nodes(:,1:3);
            nodeData.cubeCoords = zeros(size(nodeData.nodes,1),1);
            pCube = p.local(:);
            bboxBig = {pCube(:).bboxBig};
            bboxSmall = {pCube(:).bboxSmall};
            for i=1:size(nodeData.nodes,1)
                idx = find(cellfun(@(x)(all([nodeData.nodes(i,:)'>x(:,1); nodeData.nodes(i,:)'<x(:,2)])), bboxSmall));
                if isempty(idx)
                    nodeData.cubeCoords(i) = 0;
                else
                    nodeData.cubeCoords(i) = idx;
                end
            end

            % Remove all nodes outside bounding box
            toRemove = nodeData.cubeCoords == 0;
            nodeData.cubeCoords(toRemove) = [];
            nodeData.nodes(toRemove,:) = [];

            %group nodes that lie in the same cube
            nodeData.flag = zeros(size(nodeData.nodes,1),1);
            groupCount = 0;
            for i = 1:size(nodeData.nodes,1)
                if(nodeData.flag(i)==0)
                    groupCount = groupCount+1;
                    idx = nodeData.cubeCoords == nodeData.cubeCoords(i);
                    groupedNodes{groupCount}.nodes = nodeData.nodes(idx,:);
                    groupedNodes{groupCount}.cubeCoords = nodeData.cubeCoords(i);
                    nodeData.flag = nodeData.flag | idx;
                end
            end

            if exist('groupedNodes', 'var')
                % calculate local isosurfaces in global coordinates 
                for i=1:size(groupedNodes,2)
                    display(['Processing cube: ' num2str(i) '/' num2str(size(groupedNodes,2))]);
                    load(pCube(groupedNodes{i}.cubeCoords).segFile);
                    % Initalize variables 
                    segIds = zeros(1,size(groupedNodes{i}.nodes,1));
                    zeroOfCube = pCube(groupedNodes{i}.cubeCoords).bboxBig';
                    % Find ids of nodes
                    for j=1:size(groupedNodes{i}.nodes,1) 
                        rel_coords = groupedNodes{i}.nodes(j,:) - zeroOfCube(1,:) + 2;
                        segIds(j) = seg(rel_coords(1),rel_coords(2),rel_coords(3)); 
                    end
                    %delete the zeros no neuron has color black
                    segIds(segIds == 0) = [];
                    % node threshold: change threshold in last line here to apply different node threshold        
                    counts = histc(segIds,unique(segIds));
                    segIds = unique(segIds);
                    segIds(counts < 1) = [];
                    % collect supervoxel and calculate isosurfaces
                    if ~isempty(segIds)
                        cube = false(size(seg));
                        for k = 1 : size(segIds,2)
                            cube(seg == segIds(k)) = 1;
                        end
                        cube = imclose(cube, ones([3,3,3]));
                        cube = padarray(cube, [5 5 5]);
                        cube = smooth3(cube, 'gaussian', 5, 2);
                        issfs{i} = isosurface(cube, .05);
                        if ~isempty(issfs{i}.vertices)
                            issfs{i} = reducepatch(issfs{i}, .1);
                            issfs{i}.vertices(:,[1 2]) = issfs{i}.vertices(:,[2 1]); 			    
                            issfs{i}.vertices = issfs{i}.vertices + repmat(zeroOfCube(1,:) - [5 5 5],size(issfs{i}.vertices,1),1);
                            issfs{i}.vertices = issfs{i}.vertices .* repmat([11.24 11.24 28],size(issfs{i}.vertices,1),1);
                        end
                    end
                end
                if exist('issfs', 'var')
                    % Remove empty isosurfaces
                    idx = zeros(length(issfs),1);
                    for i=1:length(issfs)
                        if isempty(issfs{i}) | isempty(issfs{i}.vertices)
                            idx(i) = 1;
                        end
                    end
                    issfs(find(idx)) = [];	
                    % Save
                    exportSurfaceToAmira(issfs, [outputPath strrep(skelFile, '.nml', ['_' num2str(sk) '.issf'])]);
                    skel_data{1}.nodes(:,1:3) = skel_data{1}.nodes(:,1:3) .* repmat([11.24 11.24 28], size(skel_data{1}.nodes,1), 1);
                    %convertKnossosNmlToHoc2(skel_data, [outputPath skelFile(1:end-4)], 0, 1, 0, 0, [1 1 1]);
                end
            end
            clear groupedNodes nodeData issfs;
        end
    end

end

