function axonClasses = buildAxonSpecificityClasses(conn, axonClasses)
    % Determines specific axons for which the probability of having a
    % synapse fraction under a null model (binomial) equal or higher than
    % the observed fraction onto a specific target class is below a
    % threshold (see pThres below).
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    %% Add specificity thresholds
    % Excitatory axons
    axonClasses(1).specs = struct;
    axonClasses(1).specs.Somata.pThresh = 0.001365;
    axonClasses(1).specs.ProximalDendrite.pThresh = 0.001529;
    axonClasses(1).specs.ApicalDendrite.pThresh = 0.004332;
    axonClasses(1).specs.SmoothDendrite.pThresh = 0.000814;
    axonClasses(1).specs.AxonInitialSegment.pThresh = 0.000162;
    
    % Inhibitory axons
    axonClasses(2).specs = struct;
    axonClasses(2).specs.Somata.pThresh = 0.014222;
    axonClasses(2).specs.ProximalDendrite.pThresh = 0.011876;
    axonClasses(2).specs.ApicalDendrite.pThresh = 0.006506;
    axonClasses(2).specs.SmoothDendrite.pThresh = 0.008696;
    axonClasses(2).specs.AxonInitialSegment.pThresh = 0.000993;
    
    % Thalamocortical axons
    axonClasses(3).specs = struct;
    axonClasses(3).specs.Somata.pThresh = 0.000699;
    axonClasses(3).specs.ProximalDendrite.pThresh = 0.008995;
    axonClasses(3).specs.ApicalDendrite.pThresh = 0.001368;
    axonClasses(3).specs.SmoothDendrite.pThresh = 0.000687;
    axonClasses(3).specs.AxonInitialSegment.pThresh = 0.000486;
    
    % Corticocortical axons
    axonClasses(4).specs = struct;
    axonClasses(4).specs.Somata.pThresh = 0.001365;
    axonClasses(4).specs.ProximalDendrite.pThresh = 0.001070;
    axonClasses(4).specs.ApicalDendrite.pThresh = 0.005042;
    axonClasses(4).specs.SmoothDendrite.pThresh = 0.000940;
    axonClasses(4).specs.AxonInitialSegment.pThresh = 0.000162;
    
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
