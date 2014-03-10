function length = getPathLengthToANode( skeleton, startIdx, endIdx )

if isempty(skeleton.parameters)
    voxelsize= [1 1 1];
else
    voxelsize = [str2double(skeleton.parameters.scale.x) ...
        str2double(skeleton.parameters.scale.y) ...
        str2double(skeleton.parameters.scale.z)];
end

length = 0;
        if size(skeleton.edges,2) == 1 
            % skeleton edges seem to which dimensionality if only
            % one edge is present
            edgePos = skeleton.nodes(skeleton.edges(:,1)',1:3);
            length = length + norm((edgePos(1,:)-edgePos(2,:)).*voxelsize);
        else
            for ed=1:size(skeleton.edges,1)
                edgePos = skeleton.nodes(skeleton.edges(ed,:),1:3);
                length = length + norm((edgePos(1,:)-edgePos(2,:)).*voxelsize);
            end
        end
    end
end

end