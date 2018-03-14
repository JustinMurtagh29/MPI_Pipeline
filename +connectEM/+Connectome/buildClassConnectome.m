function [classConnectome, targetClasses] = ...
        buildClassConnectome(conn, targetClasses)
    % classConnectome, targetClasses] = ...
    %     buildClassConnectome(conn, targetClasses)
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
    
    % Get rid of synapses onto other targets. This is useful in case
    % `targetClasses` does not covert all possible target classes.
    connectome(~connectome.targetClassId, :) = [];
    
    % Build outputs
    classConnectome = accumarray( ...
        cat(2, connectome.edges(:, 1), connectome.targetClassId), ...
        connectome.synCount, [numel(conn.axons), numel(targetClasses)]);
    targetClasses = reshape(targetClasses, 1, []);
end