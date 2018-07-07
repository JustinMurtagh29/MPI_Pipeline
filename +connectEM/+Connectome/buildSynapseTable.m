function synT = buildSynapseTable(conn, synapses, varargin)
    % synT = buildSynapseTable(conn, synapses, varargin)
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.allowDuplicates = false;
    opts = Util.modifyStruct(opts, varargin{:});
    
    synT = table;
    synT.id = cell2mat(conn.connectome.synIdx);
    
    try
        synT.area = cell2mat(conn.connectomeMeta.contactArea);
    catch
        % HACK(amotta): connectEM.Connectome.prepareForSpecificityAnalysis
        % invalidates the values in `connectomeMeta` and removes this field
        % entirely to prevent use of out-dated data. That's fine, as long
        % as nobody needs information derived from it.
    end

    synT.preAggloId = repelem( ...
        conn.connectome.edges(:, 1), ...
        cellfun(@numel, conn.connectome.synIdx));

    synT.postAggloId = repelem( ...
        conn.connectome.edges(:, 2), ...
        cellfun(@numel, conn.connectome.synIdx));

    synT.isSpine = synapses.isSpineSyn(synT.id);
    
    if ~opts.allowDuplicates
        % remove synapses occuring multiple times
        % (i.e., between two or more pairs of neurites)
       [~, uniRows, uniCount] = unique(synT.id);
        synT = synT(uniRows, :);
        
        synT.occurences = accumarray(uniCount, 1);
        synT(synT.occurences > 1, :) = [];
        synT.occurences = [];
    end
end
