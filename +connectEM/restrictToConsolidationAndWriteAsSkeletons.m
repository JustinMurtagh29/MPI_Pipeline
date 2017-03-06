function restrictToConsolidationAndWriteAsSkeletons( gt, p, pT, outputFolder )

for i=1:length(pT.local)
    matchingNodes = cell(size(gt(i).segIds));
    for j=1:length(pT.local(i).trainFile)
        skelOld = skeleton(pT.local(i).trainFile{j});
        % Extract node positions and corresponding segment ids
        nodes = cellfun(@(x)x(:,1:3), skelOld.nodes, 'uni', 0);   
        nodes = cat(1, nodes{:});
        segIdsOfGTnodes = Seg.Global.getSegIds(p, nodes);
        for k=1:length(gt(i).segIds)
            temp = nodes(ismember(segIdsOfGTnodes, gt(i).segIds{k}),:);
            matchingNodes{k}(end+1:end+size(temp,1),:) = temp;
        end
    end
    matchingNodes = cellfun(@(x)unique(x, 'rows'), matchingNodes, 'uni', false);
    matchingNodes = matchingNodes(~cellfun(@isempty, matchingNodes));
    connectEM.generateSkeletonFromNodes([outputFolder 'trainingRegion' num2str(i) '.nml'], ...
        matchingNodes, strseq('Component', 1:length(matchingNodes)));
end

end
