function [conn, syn, axonClasses] = load(param, connFile, synFile)
    % [conn, syn, axonClasses] = load(param, connFile, synFile)
    %   Loads a connectome and augments it with
    %
    %   * conn.axonMeta.isThalamocortical
    %     Localgical column vector which indicates axons that were
    %     identified as being thalamocortical.
    %
    %   * conn.axonMeta.axonClass
    %     A categorical column vector which marks excitatory, inhibitory,
    %     thalamocortical and other axons.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    %% loading data
    conn = load(connFile);
    syn = load(synFile);
    
    % intersynapse distances (for detection of TC axons)
   [connDir, connName] = fileparts(connFile);
    interSynFile = sprintf('%s_intersynapse.mat', connName);
    interSynFile = fullfile(connDir, interSynFile);
    interSyn = load(interSynFile);
    
    %% mark thalamocortical axons
    conn.axonMeta.isThalamocortical = ...
        connectEM.Axon.detectThalamocorticals(conn, interSyn);
    conn.axonMeta = ...
        connectEM.Axon.completeSynapseMeta(param, conn, syn);
    
    %% label all axon classes
    axonClasses = ...
        connectEM.Connectome.buildAxonClasses(conn, 'minSynPre', 10);
    
    conn.axonMeta.axonClass(:) = {'Other'};
    conn.axonMeta.axonClass(axonClasses(2).axonIds) = {'Inhibitory'};
    conn.axonMeta.axonClass(axonClasses(1).axonIds) = {'Corticocortical'};
    conn.axonMeta.axonClass(axonClasses(3).axonIds) = {'Thalamocortical'};
    conn.axonMeta.axonClass = categorical(conn.axonMeta.axonClass);
end