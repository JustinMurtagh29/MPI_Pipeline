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
    axonClasses(1).specs.ApicalDendrite.pThresh     = 0.011567;
    axonClasses(1).specs.SmoothDendrite.pThresh     = 0.002976;
    
    % Inhibitory axons
    axonClasses(2).specs = struct;
    axonClasses(2).specs.Perisoma.pThresh           = 0.119379;
    axonClasses(2).specs.ApicalDendrite.pThresh     = 0.016054;
    axonClasses(2).specs.SmoothDendrite.pThresh     = 0.016215;
    axonClasses(2).specs.AxonInitialSegment.pThresh = 0.000823;
    
    % Thalamocortical axons
    axonClasses(3).specs = struct;
    axonClasses(3).specs.Perisoma.pThresh           = 0.017675;
    
    % Corticocortical axons
    axonClasses(4).specs = struct;
    axonClasses(4).specs.ApicalDendrite.pThresh     = 0.014586;
    axonClasses(4).specs.SmoothDendrite.pThresh     = 0.004950;
    
    %% Find specific axons
   [classConn, targetClasses] = ...
        connectEM.Connectome.buildClassConnectome(conn);
    
    for curAxonClassIdx = 1:numel(axonClasses)
        curAxonClass = axonClasses(curAxonClassIdx);
        curSpecs = curAxonClass.specs;
        
        curNullProbs = sum(classConn, 1);
        curNullProbs = curNullProbs / sum(curNullProbs);
        
        curProbs = ...
            connectEM.Specificity.calcChanceProbs( ...
                classConn, curAxonClass.axonIds, curNullProbs, ...
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
