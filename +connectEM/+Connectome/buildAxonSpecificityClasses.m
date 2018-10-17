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
    axonClasses(1).specs.ProximalDendrite.pThresh = 0.050108;
    axonClasses(1).specs.ApicalDendrite.pThresh   = 0.026136;
    axonClasses(1).specs.SmoothDendrite.pThresh   = 0.013979;
    
    % Inhibitory axons
    axonClasses(2).specs = struct;
    axonClasses(2).specs.Somata.pThresh           = 0.099645;
    axonClasses(2).specs.ProximalDendrite.pThresh = 0.115281;
    axonClasses(2).specs.ApicalDendrite.pThresh   = 0.040949;
    axonClasses(2).specs.SmoothDendrite.pThresh   = 0.010088;
    
    %% Find specific axons
   [classConn, targetClasses] = ...
        connectEM.Connectome.buildClassConnectome(conn);
    
    for curAxonClassIdx = 1:numel(axonClasses)
        curAxonClass = axonClasses(curAxonClassIdx);
        curSpecs = curAxonClass.specs;
        
        curNullProbs = ...
            connectEM.Specificity.calcFirstHitProbs( ...
                classConn(axonClasses(curAxonClassIdx).nullAxonIds, :));
        
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
