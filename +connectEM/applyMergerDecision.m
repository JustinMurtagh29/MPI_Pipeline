function skelModified = applyMergerDecision( pT )

for i=1:length(pT.local)
    for j=1:length(pT.local(i).trainFile)
        skelCompleted = skeleton(pT.local(i).trainFile{j});
        skelMergerDecision = skeleton(pT.local(i).trainFileMerger{j});
        [comments, treeIdx, nodeIdx] = skelMergerDecision.getAllComments();
        comments = strrep(comments, 'wrong', 'false');
        [~, decision, consistent] = checkForConsistency(comments);
        display('Unused comments: ');
        display(comments(~consistent));
        for k=1:length(consistent)
            % If comment can be used (e.g. consitent labels were provided)
            % If annotation is both or false remove
            % nothing needs to be done for true or duplicate
            % duplicate because the trees have already been merged
            % Remove all non consitent labels for now (e.g. lost of 2 true
            % + some false, or both mixed with false etc.)
            if (consistent(k) && (strcmp(decision(k), 'false') || strcmp(decision(k), 'both'))) ...
                    || ~consistent(k)
                    % Find node position in mergerSkel and remove in
                    % skeleton after completition
                    mergerNodePos = skelMergerDecision.nodes{treeIdx(k)}(nodeIdx(k),1:3);
                    nodeToDelete = skelCompleted.getNodesWithCoords(mergerNodePos);
                    tIdx = find(~cellfun(@isempty, nodeToDelete));
                    % Sometimes there are multiple nodes (from one or more
                    % trees in one position)
                    for l=1:length(tIdx)
                        for m=1:length(nodeToDelete{tIdx(l)})
                            skelCompleted = skelCompleted.deleteNode( ...
                                tIdx(l), nodeToDelete{tIdx(l)}(m));
                        end
                    end
            end
        end
        skelModified(i,j) = skelCompleted;
    end
end

end

function [mIds, d, consistent] = checkForConsistency(comments)
    % Checks whether label provided are consitent
    
    temp = cellfun(@(x)regexp(x, 'merger(\d\d\d).*(true|false|duplicate|both).*', 'tokens'), comments, 'uni', 0);
    assert(all(cellfun(@numel, temp) < 2));
    mIds = zeros(size(temp));
    d = cell(size(temp));
    idx = ~cellfun(@isempty, temp);
    mIds(idx) = cellfun(@(x)str2double(x{1}{1}), temp(idx));
    d(idx) = cellfun(@(x)x{1}{2}, temp(idx), 'uni', 0);
    consistent = false(size(mIds));
    for i=1:length(mIds)
        if mIds(i) ~= 0
            idx = mIds == mIds(i);
            thisD = d(idx);
            % Either: Exactly one true output, rest false or all duplicate or all both
            thisC = (sum(strcmp(thisD, 'true')) == 1 && sum(strcmp(thisD, 'false')) == length(thisD)-1) || ...
                all(strcmp(thisD, 'duplicate')) || all(strcmp(thisD, 'both'));
            if thisC
                consistent(i) = true;
            end
        end
    end
end
