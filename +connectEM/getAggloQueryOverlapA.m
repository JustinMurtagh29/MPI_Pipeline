function getAggloQueryOverlapA(p)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>

    scratchFolder = '/tmpscratch/mberning/axonQueryResults/';
    skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b'};
    skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
    skeletonFolders = [skeletonFolders {'/tmpscratch/scchr/AxonEndings/axonQueryResults/CS_MB_L4_AxonLeftQueries_nmls/'}];

    % Lookup segment ids of nodes+neighbours of nmls in all folders defined above
    [ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, skeletonFolders, false);

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

    outFile = fullfile(p.saveFolder, 'aggloState', 'axonFlightPaths.mat');
    Util.saveStruct(outFile, out);

end

