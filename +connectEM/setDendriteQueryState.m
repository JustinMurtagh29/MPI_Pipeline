function [skeletonFolders, versionSuffix, dendriteVersion,...
    dendriteVersionNew, casesToMerge] = setDendriteQueryState(state)
    % Written by
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    
    % Opportunity to define source folders and output filenames for
    % different states of axon query analysis. Can be applied on the whole
    % set of functions used during the analysis:
    % - getAggloQueryOverlapA/B
    % - getAggloQueryOverlapB_comment
    % - flightEndingOverlapRun

    
    % First run
    if strcmp(state,'1.0') 
        % source folders for flight paths
        scratchFolder = '/u/scchr/aggloState/dendQueryAnswers/';
        skeletonFolders = {'L4_dendrite_queries_27_10_2017'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '1.0';
        dendriteVersion = 'flight_01';
       
    % Commented queries of first clean run    
    elseif strcmp(state,'1.1')
        % source folders for flight paths
        scratchFolder = '/u/scchr/aggloState/dendQueryAnswers/';
        skeletonFolders = {'L4_dendrite_queries_27_10_2017'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_1.1';
        dendriteVersion = 'flight_01';
          
    % Second run
    elseif strcmp(state,'2.0') 
        % source folders for flight paths
        scratchFolder = '/u/scchr/aggloState/dendQueryAnswers/';
        skeletonFolders = {'L4_dendrite_queries_27_10_2017' 'L4_dendrite_queries_01_11_2017'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '2.0';
        dendriteVersion = 'flight_01';
       
    % Commented queries of first clean run    
    elseif strcmp(state,'2.1')
        % source folders for flight paths
        scratchFolder = '/u/scchr/aggloState/dendQueryAnswers/';
        skeletonFolders = {'L4_dendrite_queries_01_11_2017'};
        skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
        % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
        versionSuffix = '_1.1';
        dendriteVersion = 'flight_01';
      
        
    else
        error('Unknown state ''%s''', state);
    end

end

    
    