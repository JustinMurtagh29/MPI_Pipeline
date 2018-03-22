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
    synT.area = cell2mat(conn.connectomeMeta.contactArea);

    synT.preAggloId = repelem( ...
        conn.connectome.edges(:, 1), ...
        cellfun(@numel, conn.connectome.synIdx));

    synT.postAggloId = repelem( ...
        conn.connectome.edges(:, 2), ...
        cellfun(@numel, conn.connectome.synIdx));

    synT.isSpine = synapses.isSpineSyn(synT.id);
    synT.isSoma = ismember(synT.postAggloId, ...
        find(conn.denMeta.targetClass == 'Somata'));

    % remove duplicate entries
    [~, uniRows, uniCount] = unique(synT.id);
    synT = synT(uniRows, :);
    
    if ~opts.allowDuplicates
        % remove synapses occuring multiple times
        % (i.e., between two or more pairs of neurites)
        synT.occurences = accumarray(uniCount, 1);
        synT(synT.occurences > 1, :) = [];
        synT.occurences = [];
    end
end