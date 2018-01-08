function exportPathLengthSampledAgglos(param, agglos, outDir, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    opts = struct;
    opts.numAgglos = 10;
    opts.aggloMask = [];
    
    % process additional arguments
    opts = Util.modifyStruct(opts, varargin{:});
    
    if ~isempty(opts.aggloMask)
        assert(numel(agglos) == numel(opts.aggloMask));
        aggloIds = find(opts.aggloMask(:));
    else
        aggloIds = reshape(1:numel(agglos), [], 1);
    end
    
    %%
    fprintf('Calculating path lengths... ');
    pathLen = Superagglos.mstLength( ...
        agglos(aggloIds), param.raw.voxelSize);
    fprintf('done!\n');
    
    %%
    rng(0);
   [~, randIdx] = datasample( ...
        aggloIds, opts.numAgglos, ...
        'Weights', pathLen, 'Replace', false);
    
    randIds = aggloIds(randIdx);
    randLen = pathLen(randIdx);
    
    %%
    for curIdx = 1:numel(randIds)
        curId = randIds(curIdx);
        curLenUm = randLen(curIdx) / 1E3;
        
        curAgglo = agglos(curId);
        
        curSkel = skeleton();
        curSkel = curSkel.addTree( ...
            sprintf('Agglomerate %d (%.2f µm)', curId, curLenUm), ...
            curAgglo.nodes(:, 1:3), curAgglo.edges);
        curSkel = Skeleton.setParams4Pipeline(curSkel, param);
        
        curNmlName = sprintf('%d_agglomerate-%d.nml', curIdx, curId);
        curSkel.write(fullfile(outDir, curNmlName));
    end
end