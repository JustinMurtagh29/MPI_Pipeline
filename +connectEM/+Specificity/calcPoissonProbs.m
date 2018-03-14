function probs = calcPoissonProbs(classConn, axonIds, poissAxonIds)
    % probs = calcPoissonProbs(classConn, axonIds, poissAxonIds)
    %   Calculates the probability of seeing a greater or equal number of
    %   synapses onto a target class under the Poisson model.
    %
    % classConn
    %   NxM matrix with class connectome. Rows and columns correspond to
    %   axons and target classes, respectively. The entries are number of
    %   synapses.
    %
    % axonIds
    %   Optional matrix with indices of axons for which the output is
    %   calculated. Defaults to `(1:size(classConn, 1))'`.
    %
    % poissAxonIds
    %   Optional matrix with indices of axons whose synapses are considered
    %   to build the Poisson model. Defaults to `axonIds`.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % Default values
    if ~exist('axonIds', 'var') || isempty(axonIds)
        axonIds = reshape(1:size(classConn, 1), [], 1);
    end
    
    if ~exist('poissAxonIds', 'var') || isempty(poissAxonIds)
        poissAxonIds = axonIds;
    end
    
    % Poisson model probabilities
    classProbs = sum(classConn(poissAxonIds(:), :), 1);
    classProbs = classProbs ./ sum(classProbs);
    
    % Prepare output
    axonSynCount = classConn(axonIds(:), :);
    axonLambda = classProbs .* sum(axonSynCount, 2);
    
    % Calculate probabilities
    probs = 1 - arrayfun(@poisscdf, axonSynCount - 1, axonLambda);
    probs = reshape(probs, horzcat(size(axonIds), size(classConn, 2)));
end