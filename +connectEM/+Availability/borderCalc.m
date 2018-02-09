function borderCalc(param, connFile, outDir, idx)
    % Written by
    %   Kevin M. Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    blockSize = [32, 32, 16];
    info = Util.runInfo(false);
    
    %% loading data
    maxSegId = Seg.Global.getMaxSegId(param);
    seg = loadSegDataGlobal(param.seg, param.local(idx).bboxSmall);
    conn = load(connFile, 'axons', 'dendrites');
    
    %% apply axon / dendrite mappings to segmentation
    axonLookup = Agglo.buildLUT(maxSegId, conn.axons);
    dendriteLookup = Agglo.buildLUT(maxSegId, conn.dendrites);
    
    % TODO(amotta):
    % What about overlaps between axon and dendrite agglomerates? For now,
    % axons dominate over dendrites.

    % axon ids will be negative
    generalLookup = dendriteLookup;
    generalLookup(axonLookup ~= 0) = -axonLookup(axonLookup ~= 0);

    % put non-neural segment IDs in a separate class
    % non-neural segments have a positive ID → dendrites
    nonNeuralSegId = 1 + max(generalLookup(:));
    generalLookup(~generalLookup) = nonNeuralSegId;

    % apply equivalence class mapping
    seg = double(seg);
    seg(seg ~= 0) = generalLookup(seg(seg ~= 0));

    %% find edges and borders
    % TODO(amotta): Cube borders are not properly handled yet.
   [edges, ind] = connectEM.borders.codeBenedikt(seg);

    % sanity checks
    assert(size(edges, 1) == numel(ind));
    assert(issorted(edges, 2));
    
    %% split into blocks
    corrSegId = max(seg(:)) + 1;
    seg = padarray(seg, [1, 1, 1], corrSegId);
   [X, Y, Z] = ind2sub(size(seg) , ind);

    % correct for padding
    X = X - 1;
    Y = Y - 1;
    Z = Z - 1;

    blockCount = ceil(size(seg) ./ blockSize);

    blkIdx = sub2ind( ...
        blockCount, ...
        ceil(X / blockSize(1)), ...
        ceil(Y / blockSize(2)), ...
        ceil(Z / blockSize(3)));
    
    %% process blocks
    blockData = cell(blockCount);
    for curBlkIdx = 1:prod(blockCount)
        curBlkMask = (blkIdx == curBlkIdx);
        curEdges = edges(curBlkMask, :);

        if isempty(curEdges)
            blockData{curBlkIdx} = zeros(0, 4);
            continue;
        end

        % count voxels per edges
       [curEdges, ~, curEdgeId] = unique(curEdges, 'rows');
        curVxCount = accumarray(curEdgeId, 1);

        % calculate physical area
        curVxIds = ind(curBlkMask);
        curVxIds = arrayfun( ...
            @(idx) curVxIds(curEdgeId == idx), ...
            1:size(curEdges, 1), 'UniformOutput', false);

        curAreasUm = Seg.Local.physicalBorderArea2( ...
            curVxIds, curEdges, seg, param.raw.voxelSize, true);

        % set non-neural segment ID to infinity
        curEdges(curEdges(:) == nonNeuralSegId) = inf;
        blockData{curBlkIdx} = [curEdges, curVxCount, curAreasUm];
    end
    
    %% save results
    outFile = fullfile(outDir, sprintf('cube-%d.mat', idx));
    Util.save(outFile, info, blockData);
end
