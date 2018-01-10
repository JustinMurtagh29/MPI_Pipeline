function exportToAmira(param, agglos, outDir)
    % exportToAmira(param, agglos, outDir)
    %   Exports super-agglomerates to Amira for visualization. For now this
    %   is done by writing out the minimum-spanning tree of each super-
    %   agglomerate as a separate NML file. In the future this might be
    %   changed to HOC files.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    numAgglos = numel(agglos);
    numDigits = ceil(log10(numAgglos + 1));
    
    tic;
    for curIdx = 1:numAgglos
        curAgglo = agglos(curIdx);
        curNodes = curAgglo.nodes(:, 1:3);
        
        curSkel = Skeleton.fromMST(curNodes, param.raw.voxelSize);
        curSkel = Skeleton.setParams4Pipeline(curSkel, param);
        
        curOutPath = sprintf('agglomerate_%0*d', numDigits, curIdx);
        curSkel.write(fullfile(outDir, curOutPath));
        
        Util.progressBar(curIdx, numAgglos);
    end
end