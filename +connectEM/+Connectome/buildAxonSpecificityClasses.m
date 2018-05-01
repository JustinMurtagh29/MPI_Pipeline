function axonClasses = buildAxonSpecificityClasses(conn, axonClasses)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    %% Add specificity thresholds
    % Excitatory axons
    axonClasses(1).specs = struct;
    axonClasses(1).specs.ProximalDendrite.pThresh = 0.045;
    axonClasses(1).specs.ApicalDendrite.pThresh = 0.035;
    
    % Inhibitory axons
    axonClasses(2).specs = struct;
    axonClasses(2).specs.Somata.pThresh = 0.115;
    axonClasses(2).specs.ProximalDendrite.pThresh = 0.12;
    axonClasses(2).specs.ApicalDendrite.pThresh = 0.06;
    axonClasses(2).specs.SmoothDendrite.pThresh = 0.06;
    
    % Thalamocortical axons
    axonClasses(3).specs = struct;
    axonClasses(3).specs.ProximalDendrite.pThresh = 0.195;
    axonClasses(3).specs.ApicalDendrite.pThresh = 0.015;
    
    % Corticocortical axons
    axonClasses(4).specs = struct;
    axonClasses(4).specs.ProximalDendrite.pThresh = 0.03;
    axonClasses(4).specs.ApicalDendrite.pThresh = 0.04;
    axonClasses(4).specs.SmoothDendrite.pThresh = 0.04;
    
    %% Find specific axons
   [classConn, targetClasses] = ...
        connectEM.Connectome.buildClassConnectome(conn);
    
    for curAxonClassIdx = 1:numel(axonClasses)
        curAxonClass = axonClasses(curAxonClassIdx);
        curSpecs = curAxonClass.specs;
        
        curProbs = ...
            connectEM.Specificity.calcChanceProbs( ...
                classConn, ...
                curAxonClass.axonIds, ...
                curAxonClass.nullAxonIds, ...
                'distribution', 'binomial');
            
        curTargetClasses = fieldnames(curSpecs);
        for curTargetClassIdx = 1:numel(curTargetClasses)
            curTargetClass = curTargetClasses{curTargetClassIdx};
            
            curProb = curProbs(:, targetClasses == curTargetClass);
            curProbThresh = curSpecs.(curTargetClass).pThresh;
            
            curAxonIds = curProb < curProbThresh;
            curAxonIds = curAxonClass.axonIds(curAxonIds);
            curSpecs.(curTargetClass).axonIds = curAxonIds;
        end
        
        curAxonClass.specs = curSpecs;
        axonClasses(curAxonClassIdx) = curAxonClass;
    end
end
