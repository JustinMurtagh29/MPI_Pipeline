function skel = correctSkeletonsToBBox( skel, sizeCube )
% WARNING: Currently only opperates on nodesNumDataAll and nodes
for l=1:size(skel,2)
    if isfield(skel{l}, 'nodesNumDataAll') && ~isempty(skel{l}.nodesNumDataAll)
        % Correct for skeletons running out of minicube
        tooSmall = skel{l}.nodesNumDataAll(:,3:5) < ones(size(skel{l}.nodesNumDataAll,1),3);
        tooBig = skel{l}.nodesNumDataAll(:,3:5) > repmat(sizeCube,size(skel{l}.nodesNumDataAll,1),1);
        toDel = any(tooSmall,2) | any(tooBig,2);
        if sum(tooSmall(:)) || sum(tooBig(:))
            display(['Nodes to be removed: ' num2str(sum(tooSmall,1)) ' too small, ' num2str(sum(tooBig,1)) ' too big.']);
        end
        skel{l}.nodesNumDataAll(toDel,:) = [];
        edgesToDel = find(toDel);
        edgeIdx = unique(skel{l}.edges(:));
        for idx=1:length(edgesToDel)
            [row,~] = find(edgesToDel(idx) == skel{l}.edges);
	    if ~isempty(row)
		skel{l}.edges(row,:) = [];
	    end
            edgeIdx(edgeIdx == edgesToDel(idx)) = [];
        end
        edgeIdxNew = (1:length(edgeIdx))';
        for idx=1:length(edgeIdxNew)
            renumber = skel{l}.edges == edgeIdx(idx);
            skel{l}.edges(renumber) = edgeIdxNew(idx);
        end
        if sum(toDel ~= 0)
            display(['Skeleton ' num2str(l) ': ' num2str(sum(toDel)) ' nodes removed.']);
        end
    end
end

end
