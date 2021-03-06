function getDendQueryAxonAggloOverlapB(param,state)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>

    % Directory with input / output data
    dataDir = fullfile(param.saveFolder, 'aggloState');
    
    % State of query generation
    [~, suffixVersion, flighpathVersion, suffixDendrites] = connectEM.setDendriteQueryState(state);  
    
    % Load flight paths
    m = load(fullfile(dataDir, strcat('dendriteFlightPaths_',flighpathVersion,'.mat')), 'ff');
    ff = m.ff;

    % Load axon agglomerates
    m = load(fullfile(dataDir, 'axons_13_a.mat'));
    axons = m.axons(m.indBigAxons);
    axons = arrayfun(@Agglo.fromSuperAgglo, axons, 'UniformOutput', false);
    clear m;
    m = load(fullfile(dataDir, sprintf('dendrites_%s.mat',suffixDendrites)));
    dendrites = m.dendrites(m.indBigDends);
    dendrites = arrayfun(@Agglo.fromSuperAgglo, dendrites, 'UniformOutput', false);
    clear m;


    % Get flight paths without comment but with start node
    ff = structfun(@(x)x(cellfun(@isempty, ff.comments)), ff, 'uni', 0);
    ff = structfun(@(x)x(~cellfun(@isempty, ff.startNode)), ff, 'uni', 0);

    m = load(fullfile(dataDir, strcat('dendritePostQueryAnalysisState_',suffixVersion,'.mat')),'idxNoClearEnd','idxNoClearStart');
    idxNoClearEndButClearStart = find(m.idxNoClearEnd & ~m.idxNoClearStart);
    clear m
    ff = structfun(@(x)x(idxNoClearEndButClearStart), ff, 'uni', 0);
    
    segmentsLeftover = [];

    % Calculate overlap of all queries with segments
    [uniqueSegments, neighboursStartNode, ~, ~] = cellfun(@connectEM.queryAnalysis, ...
            ff.segIds, ff.neighbours, ff.nodes', ff.startNode', 'uni', 0);
    % Determine all overlaps of agglomerations with given queries
    [~, queryOverlap] = connectEM.queryAgglomerationOverlap(axons, segmentsLeftover, uniqueSegments, neighboursStartNode);
    [~, queryOverlapDendrites] = connectEM.queryAgglomerationOverlap(dendrites, segmentsLeftover, uniqueSegments, neighboursStartNode);
    % Make decision(s), here evidence/occurence threshold is applied
    % Multiple ends (all above 53vx evidence, corresponds to 2 full nodes)
    endAxon = arrayfun(@(x)x.eqClasses(x.occurences > 53), queryOverlap.ends, 'uni', 0);
    % Check dendrite attachments at start
    startDendrite = arrayfun(@(x)x.eqClasses(x.occurences > 13), queryOverlapDendrites.start, 'uni', 0);
    % Exclude all queries that do not have (at least one) clear end
    idxNoClearAxonAttachment = cellfun('isempty', endAxon);
    idxNoClearDendriteAttachment = cellfun('isempty', startDendrite);
    % 18.5% of queries excluded overall due to missing start or end (or both)
    idxGood = ~idxNoClearAxonAttachment;
    idxNoClearEndAndNoAxonAttachment = idxNoClearEndButClearStart(idxNoClearAxonAttachment);
    % Display some statistics
    display([num2str(sum(~idxNoClearAxonAttachment)./numel(idxNoClearAxonAttachment)*100, '%.2f') '% of remaining queries have axon attachment']);
    display([num2str(sum(~idxNoClearAxonAttachment)) ' in total']);
    display([num2str(sum(idxNoClearDendriteAttachment)./numel(idxNoClearDendriteAttachment)*100, '%.2f') '% of remaining queries have no start attachment at dendrite']);
    display([num2str(sum(idxNoClearDendriteAttachment)) ' in total']);
    
    results = struct;
    results.endAxon = endAxon;
    results.startDendrite = startDendrite;
    results.ff = ff;
    results.idxGood = idxGood;
    results.gitInfo = Util.gitInfo();

    % Save results and deprive writing permission
    saveFile = fullfile(dataDir, strcat('axonDendriteQueryOverlaps_',suffixVersion,'.mat'));
    save(saveFile, 'results', 'queryOverlap', 'idxNoClearAxonAttachment','idxNoClearDendriteAttachment','idxNoClearEndAndNoAxonAttachment');
    system(['chmod -w ' saveFile]);
    
    
end

