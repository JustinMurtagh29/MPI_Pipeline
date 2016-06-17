function equivMatrix = evalSegError( segmentation, skeletons)

%uniqueSegments = unique(segmentation);
%uniqueSegments(uniqueSegments == 0) = [];

equivMatrix = zeros(size(skeletons,2), single(max(segmentation(:))));
for l=1:size(skeletons,2)
    if size(skeletons{l}.nodesNumDataAll,1) > 0
        nodes = skeletons{l}.nodesNumDataAll(:,3:5);
        for m=1:size(nodes,1)
            idAtNode = segmentation(nodes(m,1), nodes(m,2), nodes(m,3));
	    %idxAtNode = find(idAtNode == uniqueSegments);
            if idAtNode
               equivMatrix(l,idAtNode) = equivMatrix(l,idAtNode) + 1;
            end
        end
    end              
end


end

