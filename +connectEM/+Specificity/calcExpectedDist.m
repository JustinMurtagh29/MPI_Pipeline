function [synFrac, expAxonCount] = ...
        calcExpectedDist(axonSynCounts, synProb, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.distribution = 'poisson';
    opts.outputFormat = 'relative';
    opts = Util.modifyStruct(opts, varargin{:});
    
    % Determine synapse number
   [synCounts, ~, synCountAxons] = unique(axonSynCounts);
    synCountAxons = accumarray(synCountAxons, 1);
    
    switch opts.distribution
        case 'poisson'
            expAxonFunc = @(nSyn, nAxons) ...
                nAxons * poisspdf((0:nSyn)', nSyn * synProb);
        case 'binomial'
            expAxonFunc = @(nSyn, nAxons) ...
                nAxons * binopdf((0:nSyn)', nSyn, synProb);
        case 'bb' % beta-binomial (marginal of dirichlet-multinomial)
            expAxonFunc = @(nSyn, nAxons) ...
                nAxons * Math.Prob.bbinopdf( ...
                    (0:nSyn)', nSyn, opts.a, opts.b);
        otherwise
            error( ...
                'Unknown probability distribution "%s"', ...
                opts.distribution);
    end
    
    t = table;
    t.expAxonCount = cell2mat(arrayfun(...
        expAxonFunc, synCounts, synCountAxons, ...
        'UniformOutput', false));
    t.synCount = cell2mat(arrayfun( ...
        @(nSyn) transpose(0:nSyn), ...
        synCounts, 'UniformOutput', false));
    t.synCountAll = reshape(repelem( ...
        synCounts(:), synCounts(:) + 1), [], 1);
    
    switch opts.outputFormat
        case 'relative'
            t.synFrac = t.synCount ./ t.synCountAll;
            t = sortrows(t, 'synFrac');
            
            % It's possible that we have multiple rows with the same
            % `synFrac`. Let's merge these entries by adding up
            % `expAxonCount`.
           [synFrac, ~, uniGroups] = unique(t.synFrac);
            expAxonCount = accumarray(uniGroups, t.expAxonCount);
        case 'absolute'
            synFrac = [t.synCount, t.synCountAll];
            expAxonCount = t.expAxonCount;
        otherwise
            error('Unknown output format "%s"', opts.outputFormat);
    end
end
