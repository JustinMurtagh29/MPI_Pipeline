function classConnectome = buildClassConnectome(conn, targetClasses)
    % classConnectome = buildClassConnectome(conn)
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % Default for `targetClasses`
    if ~exist('targetClasses', 'var') || isempty(targetClasses)
        targetClasses = unique(conn.denMeta.targetClass);
    end
    
    % Count synapses per axon
    connectome = conn.connectome;
    connectome.synCount = cellfun(@numel, connectome.synIdx);

    % Add target class to connectome
   [~, denMetaRow] = ismember(connectome.edges(:, 2), conn.denMeta.id);
    connectome.targetClass = conn.denMeta.targetClass(denMetaRow);

   [~, connectome.targetClassId] = ismember( ...
        connectome.targetClass, targetClasses);

    classConnectome = accumarray( ...
        cat(2, connectome.edges(:, 1), connectome.targetClassId), ...
        connectome.synCount, [numel(conn.axons), numel(targetClasses)]);
end