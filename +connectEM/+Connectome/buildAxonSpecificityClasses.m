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
    axonClasses(1).specs.Somata.pThresh = 0.009104;
    axonClasses(1).specs.ProximalDendrite.pThresh = 0.046053;
    axonClasses(1).specs.ApicalDendrite.pThresh = 0.071580;
    axonClasses(1).specs.SmoothDendrite.pThresh = 0.036117;
    axonClasses(1).specs.AxonInitialSegment.pThresh = 0.001264;
    
    % Inhibitory axons
    axonClasses(2).specs = struct;
    axonClasses(2).specs.Somata.pThresh = 0.145829;
    axonClasses(2).specs.ProximalDendrite.pThresh = 0.144979;
    axonClasses(2).specs.ApicalDendrite.pThresh = 0.105481;
    axonClasses(2).specs.SmoothDendrite.pThresh = 0.081366;
    axonClasses(2).specs.AxonInitialSegment.pThresh = 0.002708;
    
    % Thalamocortical axons
    axonClasses(3).specs = struct;
    axonClasses(3).specs.Somata.pThresh = 0.015878;
    axonClasses(3).specs.ProximalDendrite.pThresh = 0.151037;
    axonClasses(3).specs.ApicalDendrite.pThresh = 0.030102;
    axonClasses(3).specs.SmoothDendrite.pThresh = 0.003046;
    axonClasses(3).specs.AxonInitialSegment.pThresh = 0.019116;
    
    % Corticocortical axons
    axonClasses(4).specs = struct;
    axonClasses(4).specs.Somata.pThresh = 0.009104;
    axonClasses(4).specs.ProximalDendrite.pThresh = 0.033506;
    axonClasses(4).specs.ApicalDendrite.pThresh = 0.086021;
    axonClasses(4).specs.SmoothDendrite.pThresh = 0.047576;
    axonClasses(4).specs.AxonInitialSegment.pThresh = 0.001264;
    
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
