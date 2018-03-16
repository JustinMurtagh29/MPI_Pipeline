function [classConnectome, targetClasses, axonClasses] = ...
        buildClassConnectome(conn, varargin)
    % classConnectome, targetClasses] = ...
    %     buildClassConnectome(conn, varargin)
    %   Takes a full connectome and converts it into a class connectome by
    %   adding up the synapses within axon and / or dendrite classes.
    %     
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.axonClasses = [];
    opts.targetClasses = unique(conn.denMeta.targetClass);
    opts = Util.modifyStruct(opts, varargin{:});
    
    % Count synapses per axon
    connectome = conn.connectome;
    connectome.synCount = cellfun(@numel, connectome.synIdx);
    
    % Add axon classes to connectome
    if isempty(opts.axonClasses)
        connectome.axonClassId = connectome.edges(:, 1);
        opts.axonClasses = 1:numel(conn.axons);
    else
       [~, axonMetaRow] = ismember( ...
            connectome.edges(:, 1), conn.axonMeta.id);
        connectome.axonClassId = ...
            conn.axonMeta.axonClass(axonMetaRow);
       [~, connectome.axonClassId] = ismember( ...
            connectome.axonClassId, opts.axonClasses);
    end

    % Add target classes to connectome
    if isempty(opts.targetClasses)
        connectome.targetClassId = connectome.edges(:, 2);
        opts.targetClasses = 1:numel(conn.dendrites);
    else
       [~, denMetaRow] = ismember( ...
            connectome.edges(:, 2), conn.denMeta.id);
        connectome.targetClassId = ...
            conn.denMeta.targetClass(denMetaRow);
       [~, connectome.targetClassId] = ismember( ...
            connectome.targetClassId, opts.targetClasses);
    end
    
    % Get rid of synapses onto other targets. This is useful in case
    % `targetClasses` does not covert all possible target classes.
    % Obviously, the same holds for axons too.
    connectome(~connectome.axonClassId, :) = [];
    connectome(~connectome.targetClassId, :) = [];
    
    % Build outputs
    classConnectomeSize = ...
       [numel(opts.axonClasses), numel(opts.targetClasses)];
    classConnectome = accumarray( ...
       [connectome.axonClassId, connectome.targetClassId], ...
        connectome.synCount, classConnectomeSize);
    
    targetClasses = reshape(opts.targetClasses, 1, []);
    axonClasses = reshape(opts.axonClasses, [], 1);
end