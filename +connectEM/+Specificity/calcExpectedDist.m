function [synFrac, expAxonCount] = ...
        calcExpectedDist(axonSynCounts, synProb, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.distribution = 'poisson';
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
        otherwise
            error( ...
                'Unknown probability distribution "%s"', ...
                opts.distribution);
    end
    
    t = table;
    t.expAxonCount = cell2mat(arrayfun(...
        expAxonFunc, synCounts, synCountAxons, ...
        'UniformOutput', false));
    t.synFrac = cell2mat(arrayfun( ...
        @(nSyn) linspace(0, 1, nSyn + 1)', ...
        synCounts, 'UniformOutput', false));
    t = sortrows(t, 'synFrac');
    
    % It's possible that we have multiple rows with the same `synFrac`.
    % Let's merge these entries by adding up `expAxonCount`.
   [synFrac, ~, uniGroups] = unique(t.synFrac);
    expAxonCount = accumarray(uniGroups, t.expAxonCount);
end