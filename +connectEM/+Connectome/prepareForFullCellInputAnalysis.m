function [dendMeta, classConn] = ...
        prepareForFullCellInputAnalysis(dendMeta, classConn)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    cellT = dendMeta(dendMeta.cellId > 0, :);
    cellT(cellT.targetClass == 'AxonInitialSegment', :) = [];

   [~, cellMeta, cellT.cellId] = unique(cellT.cellId);

    cellMeta = cellT(cellMeta, :);
    cellMeta.parentId(:) = 0;
    cellMeta.targetClass(:) = {'FullInput'};
    cellMeta.synCount = accumarray(cellT.cellId, cellT.synCount);
    cellMeta.spineSynCount = accumarray(cellT.cellId, cellT.spineSynCount);
    cellMeta.id = size(classConn, 2) + transpose(1:height(cellMeta));

    cellClassConn = accumarray( ...
        cellT.cellId, cellT.id, [], @(ids) {ids});
    cellClassConn = reshape(cellClassConn, 1, []);
    cellClassConn = cell2mat(cellfun( ...
        @(ids) sum(classConn(:, ids), 2), ...
        cellClassConn, 'UniformOutput', false));

    dendMeta = cat(1, dendMeta, cellMeta);
    classConn = cat(2, classConn, cellClassConn);
end
