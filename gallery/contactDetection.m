function contactDetection(skeletons)
% Pass exactly two file names in a cell array

for id=1:2
        skel_data = readNml(skeletons{id);
        nodes = skel_data{1,1}.nodes(:,1:3);
        % for each node, find cube in which it resides
                nodeData.nodes = nodes;
                nodeData.cubeCoords = zeros(1,size(nodes,2));
                for i=1:size(nodes,1)
                   nodeData.cubeCoords(i,:) = floor(( nodes(i,:) - 1) / 128 );
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
end

end
