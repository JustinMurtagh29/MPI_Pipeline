function probs = calcChanceProbs( ...
        classConn, axonIds, nullTargetClassProbs, varargin)
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
    % nullTargetClassProbs
    %   Vector where entry `nullTargetClassProbs(i)` indicates the
    %   probability of innervating target class `i` under the null model.
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
    
    nullTargetClassProbs = reshape(nullTargetClassProbs, 1, []);
    assert(numel(nullTargetClassProbs) == size(classConn, 2));
    
    % Prepare output
    axonSynCounts = classConn(axonIds(:), :);
    
    % Calculate probabilities
    switch opts.distribution
        case 'poisson'
            probs = 1 - arrayfun( ...
                @poisscdf, axonSynCounts - 1, ...
                nullTargetClassProbs .* sum(axonSynCounts, 2));
            
        case 'binomial'
            probs = nan(size(axonSynCounts));
            for curIdx = 1:size(probs, 1)
                probs(curIdx, :) = 1 - binocdf( ...
                    axonSynCounts(curIdx, :) - 1, ...
                    sum(axonSynCounts(curIdx, :), 2), ...
                    nullTargetClassProbs);
            end
            
        case 'drmn' % dirichlet-multinomial
            % requires input 'alpha' and the drmn parameter vector
            
            assert(isfield(opts, 'alpha'), ...
                'Specify the drml alpha parameter.');
            probs = nan(size(axonSynCounts));
            % calculates the null probs separately for each target class
            % using the marginal distribution (beta-binomial)
            alpha = opts.alpha;
            n_syn_tot = sum(axonSynCounts, 2);
            for i = 1:length(alpha)
                a = alpha(i);
                b = sum(alpha(setdiff(1:length(alpha), i)));
                n_syn_class = axonSynCounts(:, i);
                for n = unique(n_syn_tot(:)')
                    % marginal is beta binomial
                    warning('off', 'all');
                    pdf = Math.Prob.bbinopdf(0:n, n, a, b);
                    warning('on', 'all');
                    ucdf = flip(cumsum(flip(pdf)));
                    thisAx = n_syn_tot == n;
                    probs(thisAx, i) = ucdf(n_syn_class(thisAx) + 1);
                end
                
            end
            assert(~any(isnan(probs(:))));
            
        otherwise
            error( ...
                'Unknown probability distribution "%s"', ...
                opts.distribution);
    end
    
    probs = reshape(probs, horzcat(size(axonIds, 1), size(classConn, 2)));
end
