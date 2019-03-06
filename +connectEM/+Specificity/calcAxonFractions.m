function [ ...
        fracSpec, ...
        fracOfSpecOntoTarget, ...
        fracOfSpecOntoTargetBars] = ...
            calcAxonFractions(axonClasses, targetClasses)
    % calcAxonFractions(axonClasses, targetClasses)
    %   Calculates the fraction of axons with target specificities.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    fracSpec = nan( ...
        numel(axonClasses), 1);
    fracOfSpecOntoTarget = nan( ...
        numel(axonClasses), numel(targetClasses));
    fracOfSpecOntoTargetBars = nan( ...
        numel(axonClasses), 2 * numel(targetClasses) - 1);

    for curId = 1:numel(axonClasses)
        curAxonClass = axonClasses(curId);
        curAxonIds = curAxonClass.axonIds;
        curSpecs = curAxonClass.specs;

        curSpecClasses = fieldnames(curSpecs);
       [~, curIds] = ismember(curSpecClasses, targetClasses);
        assert(all(curIds));

        curSpecAxonIds = cellfun( ...
            @(name) curSpecs.(name).axonIds, ...
            curSpecClasses, 'UniformOutput', false);

        % Determine fraction of axons with specificity
        curFracSpecific = mean(ismember( ...
            curAxonIds, cell2mat(curSpecAxonIds)));

        % Fraction of specific axons per target
       [curA, curB] = ndgrid(1:numel(curIds), 1:numel(curIds));
        curFracMat = zeros(numel(targetClasses));
        curFracMat(curIds, curIds) = cellfun( ...
            @(idsOne, idsTwo) numel(intersect(idsOne, idsTwo)), ...
            curSpecAxonIds(curA), curSpecAxonIds(curB));

        % NOTE(amotta): These fractions may add up to more than 100 %.
        % Reason for this is that an axon may be specific for multiple
        % target classes.
        curFracPerTarget = diag(curFracMat, 0);
        curFracPerTarget = curFracPerTarget / numel(curAxonIds);
        curFracPerTarget = curFracPerTarget / curFracSpecific;
        
        curDiag = diag(curFracMat, 0);
        curOff = diag(curFracMat, 1);

        curDiag(1:(end - 1)) = curDiag(1:(end - 1)) - curOff;
        curDiag(2:end) = curDiag(2:end) - curOff;

        curDiag = reshape(curDiag, 1, []);
        curOff = reshape(curOff, 1, []);

        curBars = cat(1, curDiag, cat(2, curOff, 0));
        curBars = curBars(1:(end - 1));
        curBars = curBars / sum(curBars);

        fracSpec(curId) = curFracSpecific;
        fracOfSpecOntoTarget(curId, :) = curFracPerTarget;
        fracOfSpecOntoTargetBars(curId, :) = curBars;
    end
end
