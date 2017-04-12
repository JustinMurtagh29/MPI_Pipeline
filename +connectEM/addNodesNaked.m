function obj = addNodesNaked(obj, tree_index, pos, diameter)
     if nargin() == 3
         diameter=1.5;
     end
     obj.nodes{tree_index} = [obj.nodes{tree_index}; pos diameter];
     % obj.nodesNumDataAll{tree_index} = [obj.nodesNumDataAll{tree_index}; obj.largestID+1 diameter pos 0 0 0 0];
     ts.id = num2str(obj.largestID+1);
     ts.radius = num2str(diameter);
     ts.x = num2str(pos(1));
     ts.y = num2str(pos(2));
     ts.z = num2str(pos(3));
     ts.inVp = '0';
     ts.inMag = '0';
     ts.time = '0';
     ts.comment = '';
     obj.nodesAsStruct{tree_index} = [obj.nodesAsStruct{tree_index} ts];
     obj.largestID = obj.largestID+1;
 end
