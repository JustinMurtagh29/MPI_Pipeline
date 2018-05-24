function skel = inputsOutputsToSkel(param, wcT, synapses, segPoints, skel)
    % skel = inputsOutputToSkel(wcT, syn, segPoints, skel)
    %   Builds a skeleton with all synaptic inputs / outputs.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    if ~exist('skel', 'var') || isempty(skel)
        skel = skeleton;
    end
    
    numDigits = ceil(log10(1 + size(wcT, 1)));
    for curIdx = 1:size(wcT, 1)
        
        curAgglo = wcT.agglo(curIdx);
        curSynIds = wcT.synapses{curIdx}.id;
        
        curPrefix = sprintf('%0*d', numDigits, curIdx);
        curName = sprintf('%s. %s', curPrefix, wcT.title{curIdx});
        
        curSynNames = arrayfun( ...
            @(id) sprintf('Synapse %d', id), ...
            curSynIds, 'UniformOutput', false);
        curSynNames = strcat( ...
            curName, {'. '}, curSynNames);
        curSynNodes = cellfun( ...
            @(a, b) segPoints(vertcat(a, b), :), ...
            synapses.presynId(curSynIds), ...
            synapses.postsynId(curSynIds), ...
            'UniformOutput', false);
        
        skel = Superagglos.toSkel(curAgglo, skel);
        skel.names{end} = curName;
        
        skel = Skeleton.fromMST(curSynNodes, param.raw.voxelSize, skel);
        skel.names((end - (numel(curSynNames) - 1)):end) = curSynNames;
    end
end
