function skel = correctSkeletonsToBBox_v2( skel, sizeCube )

skelToDel = zeros(size(skel,2),1);
overallTooSmall = zeros(1,3);
overallTooBig = zeros(1,3);
for l=1:size(skel,2)
    if isfield(skel{l}, 'nodesNumDataAll') && ~isempty(skel{l}.nodesNumDataAll)
        % Correct for skeletons running out of minicube, remove nodes
        tooSmall = skel{l}.nodesNumDataAll(:,3:5) < ones(size(skel{l}.nodesNumDataAll,1),3);
        tooBig = skel{l}.nodesNumDataAll(:,3:5) > repmat(sizeCube,size(skel{l}.nodesNumDataAll,1),1);
        toDel = any(tooSmall,2) | any(tooBig,2);
%         if sum(tooSmall(:)) || sum(tooBig(:))
%             display(['Nodes to be removed: ' num2str(sum(tooSmall,1)) ' too small, ' num2str(sum(tooBig,1)) ' too big.']);
%         end
        skel{l}.nodesNumDataAll(toDel,:) = [];
        skel{l}.nodes(toDel,:) = [];
        skel{l}.nodesAsStruct(toDel) = [];
        % ... remove edges accordingly
        if size(skel{l}.nodes,1)
            edgesToDel = find(toDel);
            edgeIdx = unique(skel{l}.edges(:));
            for idx=1:length(edgesToDel)
                [row,~] = find(edgesToDel(idx) == skel{l}.edges);
                skel{l}.edges(row,:) = [];
                edgeIdx(edgeIdx == edgesToDel(idx)) = [];
            end
            edgeIdxNew = (1:length(edgeIdx))';
            for idx=1:length(edgeIdxNew)
                renumber = skel{l}.edges == edgeIdx(idx);
                skel{l}.edges(renumber) = edgeIdxNew(idx);
            end
        else
            skelToDel(l) = 1;
            display(['Skeleton ' num2str(l) ': empty after bbox cutoff']);
        end
        if sum(toDel ~= 0)
            overallTooSmall = overallTooSmall + sum(tooSmall,1);
            overallTooBig = overallTooBig + sum(tooBig,1);
        end
    else
        skelToDel(l) = 1;
    end
end
skel = skel(~skelToDel);
display(['SUMMARY:' char(10) 'Removed ' num2str(sum(skelToDel)) ' skeletons outside bbox.' char(10) ...
    'Nodes too small (in each coordinate): ' num2str(overallTooSmall) char(10) ...
    'Nodes too big (in each coordinate): ' num2str(overallTooBig) char(10)]);

end
