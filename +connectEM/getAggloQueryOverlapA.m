function getAggloQueryOverlapA(param,state,type)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    
    % Input parameters:
    %   param
    %     Parameter structure of pipeline run
    %   state
    %     StateID of current run
    %   type
    %     Agglo type: axons=1, dendrite=0
    
    if nargin < 3
        type = 1;
    end

    % Set current state of queries
    if type ==1
        [skeletonFolders, suffix] = connectEM.setQueryState(state);    
    else
        [skeletonFolders, suffix] = connectEM.setDendriteQueryState(state);
    end

    % Lookup segment ids of nodes+neighbours of nmls in all folders defined above
    [ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments, ff.errors] = connectEM.lookupNmlMulti(param, skeletonFolders, false);

    tabulate(cellfun(@(x)x{1}{1}, cellfun(@(x)regexp(x, 'content="(.*)"', 'tokens'), ...
        cat(1, ff.comments{~cellfun(@isempty, ff.comments)}), 'uni', 0), 'uni', 0))

    % queries commented by the HIWI in terms of unsolvability (~10%)
    idx_comment = ~cellfun(@isempty,ff.comments);
    display([num2str(sum(idx_comment)) '/' num2str(numel(idx_comment)) ' queries contain comment and will not be used']);

    % 787 queries do not have a start node, not sure why (maybe the ones with more than one tree), maybe check later
    idx_startEmpty = cellfun(@isempty, ff.startNode);
    display([num2str(sum(idx_startEmpty)) '/' num2str(numel(idx_startEmpty)) ' queries do not have a start node']);

    out = struct;
    out.ff = ff;
    out.idx_comment = idx_comment;
    out.idx_startEmpty = idx_startEmpty;
    out.gitInfo = Util.gitInfo();

    % Save results and deprive writing permission
    if type ==1
        outFile = fullfile(param.saveFolder, 'aggloState', strcat('axonFlightPaths',suffix,'.mat')); 
    else
        outFile = fullfile(param.saveFolder, 'aggloState', strcat('dendriteFlightPaths',suffix,'.mat'));
    end
    Util.saveStruct(outFile, out);
    system(['chmod -w ' outFile])

end

