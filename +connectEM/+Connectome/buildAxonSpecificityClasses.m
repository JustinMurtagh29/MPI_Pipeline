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
    axonClasses(1).specs.Somata.pThresh = 0.007492;
    axonClasses(1).specs.ProximalDendrite.pThresh = 0.012747;
    axonClasses(1).specs.ApicalDendrite.pThresh = 0.011326;
    axonClasses(1).specs.SmoothDendrite.pThresh = 0.004438;
    axonClasses(1).specs.AxonInitialSegment.pThresh = 0.000363;
    
    % Inhibitory axons
    axonClasses(2).specs = struct;
    axonClasses(2).specs.Somata.pThresh = 0.033061;
    axonClasses(2).specs.ProximalDendrite.pThresh = 0.048615;
    axonClasses(2).specs.ApicalDendrite.pThresh = 0.018314;
    axonClasses(2).specs.SmoothDendrite.pThresh = 0.017851;
    axonClasses(2).specs.AxonInitialSegment.pThresh = 0.001267;
    
    % Thalamocortical axons
    axonClasses(3).specs = struct;
    axonClasses(3).specs.Somata.pThresh = 0.008984;
    axonClasses(3).specs.ProximalDendrite.pThresh = 0.040792;
    axonClasses(3).specs.ApicalDendrite.pThresh = 0.004631;
    axonClasses(3).specs.SmoothDendrite.pThresh = 0.000778;
    axonClasses(3).specs.AxonInitialSegment.pThresh = 0.002412;
    
    % Corticocortical axons
    axonClasses(4).specs = struct;
    axonClasses(4).specs.Somata.pThresh = 0.007492;
    axonClasses(4).specs.ProximalDendrite.pThresh = 0.011260;
    axonClasses(4).specs.ApicalDendrite.pThresh = 0.013997;
    axonClasses(4).specs.SmoothDendrite.pThresh = 0.004877;
    axonClasses(4).specs.AxonInitialSegment.pThresh = 0.001049;
    
    %% Find specific axons
   [classConn, targetClasses] = ...
        connectEM.Connectome.buildClassConnectome(conn);
    
    for curAxonClassIdx = 1:numel(axonClasses)
        curAxonClass = axonClasses(curAxonClassIdx);
        curSpecs = curAxonClass.specs;
        
        curNullProbs = classConn(curAxonClass.nullAxonIds, :);
        curNullProbs = curNullProbs ./ sum(curNullProbs, 2);
        curNullProbs = mean(curNullProbs, 1);
        
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
