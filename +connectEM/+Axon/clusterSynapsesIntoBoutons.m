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
    
    % HACK(amotta): For some stupid reason I had decided not to include
    % axons with no or only one synapse in the `interSyn` struct.  Now we
    % need to treat this case separately...
    synMask = cellfun(@numel, synIds) > 1;
    assert(isequal(find(synMask), interSyn.axonIds));
    
    boutonIds = cell(size(synIds));
    boutonIds(~synMask) = cellfun( ...
        @(vals) reshape(1:numel(vals), [], 1), ...
        synIds(~synMask), 'UniformOutput', false);
    
    boutonIds(synMask) = cellfun( ...
        @(a, b, c) forAxon(opts.cutoffDist, a, b, c), ...
        synIds(synMask), interSyn.synIds, interSyn.synToSynDists, ...
        'UniformOutput', false);
end

function boutonIds = forAxon(cutoff, synIds, synToSynIds, synToSyn)
    pairwiseDist = ~triu(true(size(synToSyn)));
    pairwiseDist = reshape(synToSyn(pairwiseDist), 1, []);

    links = linkage(pairwiseDist, 'average');
    boutonIds = cluster(links, 'cutoff', cutoff, 'criterion', 'distance');
    
    % Fix order
   [~, fixIds] = ismember(synIds, synToSynIds);
    boutonIds = reshape(boutonIds(fixIds), [], 1);
end
