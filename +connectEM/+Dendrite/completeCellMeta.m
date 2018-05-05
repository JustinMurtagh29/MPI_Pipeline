function dendMeta = completeCellMeta(param, conn)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    dendMeta = conn.denMeta;
    
    %% Add information about whole cells
    dend = conn.info.param.dendriteFile;
    dend = load(dend, 'idxWholeCells', 'idxSomata', 'idxAIS');
    
    dendMeta.cellId = horzcat( ...
        dend.idxWholeCells(dendMeta.parentId), ...
        dend.idxSomata(dendMeta.parentId), ...
        dend.idxAIS(dendMeta.parentId));
    dendMeta.cellId = max(dendMeta.cellId, [], 2);
    
    %% Add flag if part of interneuron
    segPoints = Seg.Global.getSegToPointMap(param);
    voxelSize = param.raw.voxelSize;
    
    soma = dendMeta(dendMeta.targetClass == 'Somata', :);
    soma.pos = cell2mat(cellfun( ...
        @(segIds) mean(segPoints(segIds, :), 1), ...
        conn.dendrites(soma.id), 'UniformOutput', false));
    
    % The interneuron positions were extracted from the KAMIN list:
    % https://docs.google.com/spreadsheets/d/1DVk6wqlg6bI0XTToWfMVGJ0W7XFnNmGLs3XJe4cOKtY
    inSomaPos = 1 + [ ...
        447, 4609, 2507;  % KAMIN cell 23
        3872, 5306, 395]; % KAMIN cell 31
   [~, inCellIds] = pdist2( ...
       voxelSize .* soma.pos, ...
       voxelSize .* inSomaPos, ...
       'squaredeuclidean', 'Smallest', 1);
    inCellIds = soma.cellId(inCellIds);
    
    dendMeta.isInterneuron = ...
        dendMeta.targetClass == 'SmoothDendrite' ...
      | ismember(dendMeta.cellId, inCellIds);
end
