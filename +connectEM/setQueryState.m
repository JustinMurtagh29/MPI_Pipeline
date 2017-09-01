function [skeletonFolders, fileFlightPathsSuffix, fileQueryOverlapsSuffix] = setQueryState(state)
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
        scratchFolder = '/tmpscratch/mberning/axonQueryResults/';
        skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        skeletonFolders = [skeletonFolders {'/tmpscratch/scchr/AxonEndings/axonQueryResults/CS_MB_L4_AxonLeftQueries_nmls/'}];
        % filename for flight paths in getAggloQueryOverlapA
        fileFlightPathsSuffix = '';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        fileQueryOverlapsSuffix = '';
       
    % Commented queries of first clean run    
    elseif strcmp(state,'1.1')
        % source folders for flight paths
        skeletonFolders = [];
        % filename for flight paths in getAggloQueryOverlapA
        fileFlightPathsSuffix = '';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        fileQueryOverlapsSuffix = '_1.1';
          
    % Second run
    elseif strcmp(state,'2.0')
        % source folders for flight paths
        scratchFolder = '/tmpscratch/mberning/axonQueryResults/';
        skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        skeletonFolders = [skeletonFolders {'/tmpscratch/scchr/AxonEndings/axonQueryResults/CS_MB_L4_AxonLeftQueries_nmls/'}];
        skeletonFolders = [skeletonFolders {'/tmpscratch/.../CS_MB_L4_axonEndingQueries_30_08_2017'}];
        % filename for flight paths in getAggloQueryOverlapA
        fileFlightPathsSuffix = '_2.0';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        fileQueryOverlapsSuffix = '_2.0';  
     
    % Commented queries of second run
    elseif strcmp(state,'2.1')
        % source folders for flight paths
        skeletonFolders = {'/tmpscratch/.../CS_MB_L4_axonEndingQueries_30_08_2017'};
        % filename for flight paths in getAggloQueryOverlapA
        fileFlightPathsSuffix = '_2.1';
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        fileQueryOverlapsSuffix = '_2.1';  

    end
    
    