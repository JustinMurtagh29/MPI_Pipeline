function getAggloQueryOverlapB_comment(param,state)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>

    % Directory with input / output data
    dataDir = fullfile(param.saveFolder, 'aggloState');
    
    % State of query generation
    [skeletonFolders, suffixFlightPaths, suffixOverlaps] = connectEM.setQueryState(state);   

    % Load flight paths
    m = load(fullfile(dataDir, strcat('axonFlightPaths',suffixFlightPaths,'.mat')), 'ff');
    ff = m.ff;

    % Load axon agglomerates
    m = load(fullfile(dataDir, 'axons_04.mat'));
    axons = m.axons(m.indBigAxons);
    axons = arrayfun(@Agglo.fromSuperAgglo, axons, 'UniformOutput', false);
    clear m;

    % Get flight paths with comment
    ff = structfun(@(x)x(~cellfun(@isempty, ff.comments)), ff, 'uni', 0);
    ff = structfun(@(x)x(~cellfun(@isempty, ff.startNode)), ff, 'uni', 0);

    segmentsLeftover = [];

    % Calculate overlap of all queries with segments
    [uniqueSegments, neighboursStartNode, nodesExcludedIdx, startNodeIdx] = cellfun(@connectEM.queryAnalysis, ...
            ff.segIds, ff.neighbours, ff.nodes', ff.startNode', 'uni', 0);
    % Determine all overlaps of agglomerations with given queries
    [partition, queryOverlap] = connectEM.queryAgglomerationOverlap(axons, segmentsLeftover, uniqueSegments, neighboursStartNode);
    % Make decision(s), here evidence/occurence threshold is applied
    % Always one (or none if evidence below 14, 1/2 node) start eqClass
    startAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 13), queryOverlap.start, 'uni', 0);
    % Exclude all queries that do not have a clear starting point
    idxNoClearStart = cellfun('isempty', startAgglo);
    % Multiple ends (all above 53vx evidence, corresponds to 2 full nodes)
    endAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 53), queryOverlap.ends, 'uni', 0);
    % Exclude startAgglo from endAgglo (as we do not want to count self-attachment)
    endAgglo = cellfun(@(x,y)setdiff(x,y), endAgglo, startAgglo, 'uni', 0);
    % Exclude all queries that do not have (at least one) clear end
    idxNoClearEnd = cellfun('isempty', endAgglo);
    % 18.5% of queries excluded overall due to missing start or end (or both)
    idxGood = ~(idxNoClearStart | idxNoClearEnd);
    % Display some statistics
    display([num2str(sum(idxNoClearStart)./numel(idxNoClearStart)*100, '%.2f') '% of remaining queries have no clear start']);
    display([num2str(sum(idxNoClearStart)) ' in total']);
    display([num2str(sum(idxNoClearEnd)./numel(idxNoClearEnd)*100, '%.2f') '% of remaining queries have no clear end']);
    display([num2str(sum(idxNoClearEnd)) ' in total']);
    display([num2str(sum(idxGood)./numel(idxGood)*100, '%.2f') '% of remaining queries have clear start and ending']);
    display([num2str(sum(idxGood)) ' in total']);
    display([num2str(numel(cat(2, endAgglo{idxGood}))) ' attachments made by ' num2str(sum(idxGood)) ' queries']);
    % Find CC of eqClasses to be joined including single eqClasses with or
    % without dangling query
    edgesCC = cellfun(@(x,y)combnk([x y], 2), startAgglo(idxGood), endAgglo(idxGood), 'uni', 0);
    edgesCC = cat(1,edgesCC{:});
    edgesCC = sort(edgesCC, 2);
    eqClassCC = Graph.findConnectedComponents(edgesCC, true, true);
    sizeEqClassCC = sort(cellfun(@numel, eqClassCC), 'descend');
    eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(axons), cell2mat(eqClassCC)))'];
    display(sizeEqClassCC(1:10));

    results = struct;
    results.startAgglo = startAgglo;
    results.endAgglo = endAgglo;
    results.ff = ff;
    results.idxGood = idxGood;
    results.gitInfo = Util.gitInfo();

    % Save and deprive writing permission
    saveFile = fullfile(dataDir, strcat('axonPostQueryAnalysisState',suffixOverlaps,'.mat'));
    save(saveFile);
    system(['chmod -w ' saveFile])

end

