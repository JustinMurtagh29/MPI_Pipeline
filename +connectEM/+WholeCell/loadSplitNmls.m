function nmlT = loadSplitNmls(splitNmlDir)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    nmlT = table;
    nmlT.path = dir(fullfile(splitNmlDir, '*.nml'));
    nmlT.path = reshape({nmlT.path.name}, [], 1);
    nmlT.path = fullfile(splitNmlDir, nmlT.path);
    
    nmlT.nml = cellfun(@slurpNml, nmlT.path);
    
    nmlT.aggloId = arrayfun(@(nml) regexp( ...
            nml.parameters.experiment.description, ...
            'Agglomerate (\d+)$', 'tokens', 'once'), nmlT.nml);
    nmlT.aggloId = str2double(nmlT.aggloId);

    nmlT.dendNodes = cell(size(nmlT.path));
    for curIdx = 1:numel(nmlT.dendNodes)
        curNml = nmlT.nml(curIdx);
        
        curTrees = NML.buildTreeTable(curNml);
        curNodes = NML.buildNodeTable(curNml);
        curComments = NML.buildCommentTable(curNml);

        curNodes.nodeId(:) = {''};
       [curMask, curIds] = ismember(curNodes.id, curComments.node);
        curNodes.nodeId(curMask) = curComments.comment(curIds(curMask));

        curNodes.nodeId = cellfun( ...
            @(nodeId) regexp( ...
                nodeId, 'Node (\d+)$', 'tokens', 'once'), ...
            curNodes.nodeId, 'UniformOutput', false);
        curNodes.nodeId(cellfun(@isempty, curNodes.nodeId)) = {{'0'}};
        curNodes.nodeId = str2double(vertcat(curNodes.nodeId{:}));
        curNodes(~curNodes.nodeId, :) = [];

        % Remove soma tree(s)
        % NOTE(amotta): Feel free to remove
        curSomaTreeIds = curTrees.id(contains( ...
            curTrees.name, 'soma', 'IgnoreCase', true));
        curNodes(ismember(curNodes.treeId, curSomaTreeIds), :) = [];

       [~, ~, curNodes.dendId] = unique(curNodes.treeId);
        nmlT.dendNodes{curIdx} = accumarray( ...
            curNodes.dendId, curNodes.nodeId, ...
            [max([1; curNodes.dendId(:)]), 1], ...
            @(ids) {ids}, {zeros(0, 1)});
    end
end
