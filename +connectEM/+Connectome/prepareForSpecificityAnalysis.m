function conn = prepareForSpecificityAnalysis(conn)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    %% Configuration
    targetClassOrder = { ...
        'Somata', ...
        'ProximalDendrite', ...
        'SmoothDendrite', ...
        'ApicalDendrite', ...
        'AxonInitialSegment', ...
        'OtherDendrite'};
    
    %% Treat soma-based reconstructions
    inMask = conn.denMeta.isInterneuron;
    wcMask = conn.denMeta.targetClass == 'WholeCell';

    % Inhibitory whole cell → smooth dendrite
    conn.denMeta.targetClass(wcMask &  inMask) = 'SmoothDendrite';
    
    % Excitatory whole cell → proximal dendrite
    conn.denMeta.targetClass(wcMask & ~inMask) = 'ProximalDendrite';
    
    %% Fix order of categories
    cats = unique(conn.denMeta.targetClass);
    conn.denMeta.targetClass = categorical( ...
        conn.denMeta.targetClass, cats);
    conn.denMeta.targetClass = reordercats( ...
        conn.denMeta.targetClass, targetClassOrder);
end
