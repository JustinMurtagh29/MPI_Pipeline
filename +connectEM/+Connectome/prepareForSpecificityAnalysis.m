function conn = prepareForSpecificityAnalysis(conn)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    inMask = conn.denMeta.isInterneuron;
    wcMask = conn.denMeta.targetClass == 'WholeCell';

    % Inhibitory whole cell → smooth dendrite
    conn.denMeta.targetClass(wcMask &  inMask) = 'SmoothDendrite';
    
    % Excitatory whole cell → proximal dendrite
    conn.denMeta.targetClass(wcMask & ~inMask) = 'ProximalDendrite';
end
