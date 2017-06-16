function obj = intermediateNodes(obj, tree_index, max_dist)
     % interpolate trees (inject intermediate nodes until all edge
     % lengths > max_dist
     numedges = size(obj.edges{tree_index},1);
     for i = 1:numedges
         thisEdge = obj.edges{tree_index}(i,:);
         p1 = obj.nodes{tree_index}(thisEdge(1),1:3);
         p2 = obj.nodes{tree_index}(thisEdge(2),1:3);
         dist = sqrt(sum(((p1-p2).*obj.scale).^2));
         if dist>max_dist
             numNodesToAdd = ceil(dist/max_dist)-1;
             obj.edges{tree_index}(i,2) = length(obj.nodesAsStruct{tree_index})+1;
             for j = 1:numNodesToAdd
                 nodeTodo = p1+(p2-p1)*(j-1)/numNodesToAdd;
                 obj = connectEM.addNodesNaked(obj,tree_index,round(nodeTodo));
                 if j>1
                     obj.edges{tree_index}(end+1,:) = [length(obj.nodesAsStruct{tree_index})-1 length(obj.nodesAsStruct{tree_index})];
                 end
             end
             obj.edges{tree_index}(end+1,:) = [length(obj.nodesAsStruct{tree_index}),thisEdge(2)];
         end
     end
 end
