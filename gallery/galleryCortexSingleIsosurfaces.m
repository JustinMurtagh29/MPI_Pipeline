function galleryCortexSingleIsosurfaces(p, skelPath, skelFile, outputPath)

    % Display status
    display(['Processing skeleton: ' skelFile]);
    % evalc to supress output of parseNml
    [~,skel_data] = evalc('parseNml_noInVP([skelPath skelFile])'); 

    % Find X,Y,Z for each node if it is within parameter.local(X,Y,Z).bboxSmall
    for sk = 1:length(skel_data)
        idxIso = 1;
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

            % calculate local isosurfaces in global coordinates 
            for i=1:size(groupedNodes,2)
                display(['Processing cube: ' num2str(i) '/' num2str(size(groupedNodes,2))]);
                load(pCube(groupedNodes{i}.cubeCoords).segFile);
                % Initalize variables 
                segIds = zeros(1,size(groupedNodes{i}.nodes,1));
                zeroOfCube = pCube(groupedNodes{i}.cubeCoords).bboxBig';
                % Find ids of nodes
                for j=1:size(groupedNodes{i}.nodes,1) 
                    rel_coords = groupedNodes{i}.nodes(j,:) - zeroOfCube(1,:);
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
                    for k = 1 : size(segIds,2)
                        cube = false(size(seg));
                        cube(seg == segIds(k)) = 1;
                        cube = imclose(cube, ones([3,3,3]));
                        cube = padarray(cube, [2 2 2]);
                        cube = smooth3(cube, 'gaussian', 5, 2);
                        issfs{idxIso} = isosurface(cube, .1);
                        if ~isempty(issfs{idxIso}.vertices)
                            issfs{idxIso} = reducepatch(issfs{idxIso}, .01);
                            issfs{idxIso}.vertices(:,[1 2]) = issfs{idxIso}.vertices(:,[2 1]); 			    
                            issfs{idxIso}.vertices = issfs{idxIso}.vertices + repmat(zeroOfCube(1,:) - [2 2 2],size(issfs{idxIso}.vertices,1),1);
                            issfs{idxIso}.vertices = issfs{idxIso}.vertices .* repmat([11.24 11.24 28],size(issfs{idxIso}.vertices,1),1);
                        end
                        idxIso = idxIso + 1;
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
            end
            clear groupedNodes nodeData issfs;
        end
    end

end

