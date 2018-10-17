function boutonIds = clusterSynapsesIntoBoutons(synIds, interSyn, varargin)
    % boutonIds = clusterSynapsesIntoBoutons(synIds, interSyn, varargin)
    %   Clusters the synapses in `synIds` into boutons and returns in
    %   `boutonIds{i}(j)` the bouton to which synapse `j` of axon `i`
    %   belongs.
    %
    % See also
    %   connectEM.Axon.getSynapses
    %   connectEM.Axon.Script.calculateSynToSynDists
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    
    % NOTE(amotta): For a derivation of this default value, see
    % connectEM.Axon.Script.detectMultiHitBoutons.
    opts.cutoffDist = 2433;
    opts = Util.modifyStruct(opts, varargin{:});
    
    boutonIds = cellfun( ...
        @(a, b, c) forAxon(opts.cutoffDist, a, b, c), ...
        synIds, interSyn.synIds, interSyn.synToSynDists, ...
        'UniformOutput', false);
end

function boutonIds = forAxon(cutoff, synIds, synToSynIds, synToSyn)
    if numel(synIds) < 2
        boutonIds = 1:numel(synIds);
        boutonIds = reshape(boutonIds, [], 1);
        return;
    end
    
    pairwiseDist = ~triu(true(size(synToSyn)));
    pairwiseDist = reshape(synToSyn(pairwiseDist), 1, []);
    
    links = linkage(pairwiseDist, 'average');
    boutonIds = cluster(links, 'cutoff', cutoff, 'criterion', 'distance');
    
    % Fix order
   [~, fixIds] = ismember(synIds, synToSynIds);
    boutonIds = reshape(boutonIds(fixIds), [], 1);
end
