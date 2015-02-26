function galleryCortexEmptyNodes(p, skelPath, skelFile, outputPath)

    % Display status
    display(['Processing skeleton: ' skelFile]);
    % evalc to supress output of parseNml
    [~,skel_data] = evalc('parseNml([skelPath skelFile])'); 

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
                emptyNodes(i,1) = length(segIds);
                %delete the zeros no neuron has color black
                emptyNodesDistribution{i} = segIds == 0;
                segIds(segIds == 0) = [];
                emptyNodes(i,2) = length(segIds);
            end
            % Save
            save([outputPath strrep(skelFile, '.nml', ['_' num2str(sk) '.mat'])], 'emptyNodes', 'emptyNodesDistribution');
        end
        clear groupedNodes nodeData;
    end
end

