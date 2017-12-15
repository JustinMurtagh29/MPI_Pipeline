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

    switch state
        % First run
        case '1.0'
            % source folders for flight paths
            scratchFolder = '/u/scchr/aggloState/dendQueryAnswers/';
            skeletonFolders = {'L4_dendrite_queries_27_10_2017'};
            skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
            % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
            versionSuffix = '1.0';
            dendriteVersion = 'flight_01';
            
            % Commented queries of first clean run
        case '1.1'
            % source folders for flight paths
            scratchFolder = '/u/scchr/aggloState/dendQueryAnswers/';
            skeletonFolders = {'L4_dendrite_queries_27_10_2017'};
            skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
            % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
            versionSuffix = '_1.1';
            dendriteVersion = 'flight_01';
            
            % Second run
        case '2.0'
            % source folders for flight paths
            scratchFolder = '/u/scchr/aggloState/dendQueryAnswers/';
            skeletonFolders = {'L4_dendrite_queries_27_10_2017' 'L4_dendrite_queries_01_11_2017'};
            skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
            % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
            versionSuffix = '2.0';
            dendriteVersion = 'flight_01';
            dendriteVersionNew = 'flight_02';
            
            % Commented queries of first clean run
        case '2.1'
            % source folders for flight paths
            scratchFolder = '/u/scchr/aggloState/dendQueryAnswers/';
            skeletonFolders = {'L4_dendrite_queries_01_11_2017'};
            skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
            % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
            versionSuffix = '_2.1';
            dendriteVersion = 'flight_01';
            
            % Run with the whole cells included
        case '2.2'
            % source folders for flight paths
            scratchFolder = '/u/scchr/aggloState/dendQueryAnswers/';
            skeletonFolders = {'L4_dendrite_queries_27_10_2017' 'L4_dendrite_queries_01_11_2017'};
            skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
            % filename additionals for getAggloQueryOverlapB and flightEndingOverlapRun
            versionSuffix = '2.2';
            dendriteVersion = 'andWholeCells_01';
            dendriteVersionNew = 'andWholeCells_02';
        otherwise
            error('Unknown state ''%s''', state);
    end

end

    
    