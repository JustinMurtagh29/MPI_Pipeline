function [synFrac, expAxonCount] = calcPoissonDist(axonSynCounts, synProb)
    % Determine synapse number
   [synCounts, ~, synCountAxons] = unique(axonSynCounts);
    synCountAxons = accumarray(synCountAxons, 1);

    % Poisson
    poiss = table;
    poiss.expAxonCount = cell2mat(arrayfun( ...
        @(nSyn, nAxons) nAxons * poisspdf((0:nSyn)', nSyn * synProb), ...
        synCounts, synCountAxons, 'UniformOutput', false));
    poiss.synFrac = cell2mat(arrayfun( ...
        @(nSyn) (0:nSyn)' ./ nSyn, synCounts, 'UniformOutput', false));
    poiss = sortrows(poiss, 'synFrac');
    
    % Prepare output
    synFrac = poiss.synFrac;
    expAxonCount = poiss.expAxonCount;
end