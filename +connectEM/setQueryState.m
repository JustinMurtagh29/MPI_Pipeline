function [skeletonFolders, flightPathsSuffix, versionSuffix, axonVersion,...
    axonVersionNew, casesToMerge] = setQueryState(state)
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
        axonVersion = [];
     
    % Commented queries of second run
    elseif strcmp(state,'2.1')
        % source folders for flight paths
        skeletonFolders = {'/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/CS_MB_L4_axonEndingQueries_30_08_2017_nmls/'};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '_2.1';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_2.1'; 
        axonVersion = [];

            % Third run
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
        axonVersion = [];
     
    % Commented queries of second run
    elseif strcmp(state,'3.1')
        % source folders for flight paths
        skeletonFolders = {'/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/CS_MB_L4_axonEndingQueries_30_08_2017/'};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '_3.1';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_3.1'; 
        
        % Fourth run
    elseif strcmp(state,'4.0')
        % source folders for flight paths
        scratchFolder = '/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/';
        skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b' ...
            'CS_MB_L4_AxonLeftQueries_nmls' 'CS_MB_L4_axonEndingQueries_30_08_2017' ...
            'CS_MB_L4_axEndQuerySpecial_14_09_2017'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '_4.0';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_4.0';
        axonVersion = 4;
     
    % Commented queries of second run
    elseif strcmp(state,'4.1')
        % source folders for flight paths
        skeletonFolders = {'/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/CS_MB_L4_axEndQuerySpecial_14_09_2017/'};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '_4.1';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_4.1'; 
        
     % Fourth run
    elseif strcmp(state,'5.0')
        % source folders for flight paths
        scratchFolder = '/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/';
        skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b' ...
            'CS_MB_L4_AxonLeftQueries_nmls' 'CS_MB_L4_axonEndingQueries_30_08_2017' ...
            'CS_MB_L4_axEndQuerySpecial_14_09_2017' 'CS_MB_L4_axEndQuerySpecial2_16_09_2017'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '_5.0';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_5.0';
        axonVersion = 4;
        axonVersionNew = '05';
        
    % Commented queries of second run
    elseif strcmp(state,'5.1')
        % source folders for flight paths
        skeletonFolders = {'/u/mberning/results/pipeline/20170217_ROI/aggloState/queryAnswers/CS_MB_L4_axEndQuerySpecial2_16_09_2017/'};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = '_5.1';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_5.1'; 
        
    %% Choose cases for different axon states
    elseif strcmp(state,'5.2')
        % source folders for flight paths
        skeletonFolders = {};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = [];
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_5.0'; 
        axonVersion = 4;
        axonVersionNew = '05_a';
        casesToMerge = [1:6,8];
    
    elseif strcmp(state,'5.3')
        % source folders for flight paths
        skeletonFolders = {};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = [];
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_5.0'; 
        axonVersion = 4;
        axonVersionNew = '05_b';
        casesToMerge = [1:4,6,8];
    
    elseif strcmp(state,'5.4')
        % source folders for flight paths
        skeletonFolders = {};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = [];
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_5.0'; 
        axonVersion = 4;
        axonVersionNew = '05_c';
        casesToMerge = [1:4,6,8:14];
        
    elseif strcmp(state,'5.5')
        % source folders for flight paths
        skeletonFolders = {};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = [];
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_5.0'; 
        axonVersion = 4;
        axonVersionNew = '05_d';
        casesToMerge = [1:6,8:14];
     
    elseif strcmp(state,'5.6')
        % source folders for flight paths
        skeletonFolders = {};
        % filename for flight paths in getAggloQueryOverlapA
        flightPathsSuffix = [];
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_5.0'; 
        axonVersion = 4;
        axonVersionNew = '05_E3a_';
        casesToMerge = [1:6,8:14];
        
    end

end

    
    