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
    
    try
        connSynFile = conn.info.param.synFile;
    catch
        % NOTE(amotta): `synFile` was not stored in connectome. So, let's
        % use the user specified one. We expect the user to have specified
        % a `synFile`.
        connSynFile = synFile;
    end
    
    if ~exist('synFile', 'var') || isempty(synFile)
        synFile = connSynFile;
    elseif ~strcmp(synFile, connSynFile)
        % Show an error if the user-specified synapse file doesn't match
        % the one used to generate the connectome.
        error('Specified synapse file does not match connectome');
    end
    
    syn = load(synFile);
    
    % intersynapse distances (for detection of TC axons)
   [connDir, connName] = fileparts(connFile);
    interSynFile = sprintf('%s_intersynapse_v2.mat', connName);
    interSynFile = fullfile(connDir, interSynFile);
    interSyn = load(interSynFile);
    
    % meta data about axonal boutons (for detection of TC axons)
    boutonMetaFile = sprintf('%s_axonalBoutons_v1.mat', connName);
    boutonMetaFile = fullfile(connDir, boutonMetaFile);
    boutonMeta = load(boutonMetaFile);
    
    %% complete axon meta data
    conn.axonMeta = ...
        connectEM.Axon.completeSynapseMeta( ...
            param, interSyn, boutonMeta, conn, syn);
    
    %% label dendrites
    conn.denMeta = connectEM.Dendrite.completeCellMeta(param, conn);
    
    %% label all axon classes
   [axonClasses, conn] = ...
        connectEM.Connectome.buildAxonClasses(conn, 'minSynPre', 10);
end
