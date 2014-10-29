function equivMatrixBinary = evalSegError( segmentation, skeletons, nodeThres )

nrObjects = double(max(segmentation(:)));
equivMatrix = sparse(size(skeletons,2), nrObjects);
for l=1:size(skeletons,2)
	nodes = skeletons{l}.nodes(:,1:3);
        for m=1:size(nodes,1)
		if segmentation(nodes(m,1), nodes(m,2), nodes(m,3))	
			equivMatrix(l,segmentation(nodes(m,1), nodes(m,2), nodes(m,3))) = ...
				equivMatrix(l,segmentation(nodes(m,1), nodes(m,2), nodes(m,3))) + 1;
		end
	end
end
equivMatrixBinary = equivMatrix >= nodeThres;

end

