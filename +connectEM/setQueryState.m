function [skeletonFolders, flightPathsSuffix, versionSuffix, axonVersion] = setQueryState(state)
    % Written by
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    
    % Opportunity to define source folders and output filenames for
    % different states of axon query analysis. Can be applied on the whole
    % set of functions used during the analysis:
    % - getAggloQueryOverlapA/B
    % - getAggloQueryOverlapB_comment
    % - flightEndingOverlapRun

    
    % First 'clean' run
    if strcmp(state,'1.0') 
        % source folders for flight paths
        scratchFolder = '/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/';
        skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b' ...
            'CS_MB_L4_AxonLeftQueries_nmls'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '';
        axonVersion = [];
       
    % Commented queries of first clean run    
    elseif strcmp(state,'1.1')
        % source folders for flight paths
        skeletonFolders = [];
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_1.1';
        axonVersion = [];
          
    % Second run
    elseif strcmp(state,'2.0')
        % source folders for flight paths
        scratchFolder = '/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/';
        skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b' ...
            'CS_MB_L4_AxonLeftQueries_nmls' 'CS_MB_L4_axonEndingQueries_30_08_2017_nmls'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '_2.0';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_2.0';
        axonVersion = 4;
     
    % Commented queries of second run
    elseif strcmp(state,'2.1')
        % source folders for flight paths
        skeletonFolders = {'/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/CS_MB_L4_axonEndingQueries_30_08_2017_nmls/'};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '_2.1';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_2.1'; 
        axonVersion = [];

            % Second run
    elseif strcmp(state,'3.0')
        % source folders for flight paths
        scratchFolder = '/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/';
        skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b' ...
            'CS_MB_L4_AxonLeftQueries_nmls' 'CS_MB_L4_axonEndingQueries_30_08_2017'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '_3.0';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_3.0';
        axonVersion = 4;
     
    % Commented queries of second run
    elseif strcmp(state,'3.1')
        % source folders for flight paths
        skeletonFolders = {'/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/CS_MB_L4_axonEndingQueries_30_08_2017/'};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '_3.1';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_3.1'; 
        axonVersion = [];

    end
    
    