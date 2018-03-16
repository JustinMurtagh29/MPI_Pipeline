function probs = calcChanceProbs( ...
        classConn, axonIds, nullAxonIds, varargin)
    % probs = calcChanceProbs(classConn, axonIds, poissAxonIds, varargin)
    %   Calculates the probability of seeing a greater or equal number of
    %   synapses onto a target class under a null hypothesis.
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
    % nullAxonIds
    %   Optional matrix with indices of axons whose synapses are considered
    %   to build the null model. Defaults to `axonIds`.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.distribution = 'poisson';
    opts = Util.modifyStruct(opts, varargin{:});
    
    % Default values
    if ~exist('axonIds', 'var') || isempty(axonIds)
        axonIds = reshape(1:size(classConn, 1), [], 1);
    end
    
    if ~exist('nullAxonIds', 'var') || isempty(nullAxonIds)
        nullAxonIds = axonIds;
    end
    
    % Probabilities for null model
    classProbs = sum(classConn(nullAxonIds(:), :), 1);
    classProbs = classProbs ./ sum(classProbs);
    
    % Prepare output
    axonSynCounts = classConn(axonIds(:), :);
    
    % Calculate probabilities
    switch opts.distribution
        case 'poisson'
            probs = 1 - arrayfun( ...
                @poisscdf, axonSynCounts - 1, ...
                classProbs .* sum(axonSynCounts, 2));
            
        case 'binomial'
            probs = nan(size(axonSynCounts));
            for curIdx = 1:size(probs, 1)
                probs(curIdx, :) = 1 - binocdf( ...
                    axonSynCounts(curIdx, :) - 1, ...
                    sum(axonSynCounts(curIdx, :), 2), ...
                    classProbs);
            end
            
        otherwise
            error( ...
                'Unknown probability distribution "%s"', ...
                opts.distribution);
    end
    
    probs = reshape(probs, horzcat(size(axonIds), size(classConn, 2)));
end