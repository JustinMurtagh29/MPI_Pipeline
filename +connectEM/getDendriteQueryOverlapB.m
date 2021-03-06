function getDendriteQueryOverlapB(param,state)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>

    % Directory with input / output data
    dataDir = fullfile(param.saveFolder, 'aggloState');
    segmentMeta = load(fullfile(param.saveFolder, 'segmentMeta.mat'),'voxelCount');
    % State of query generation
    [~, suffixVersion, flighpathVersion, suffixDendrites] = connectEM.setDendriteQueryState(state);    

    % Load flight paths
    m = load(fullfile(dataDir, strcat('dendriteFlightPaths_',flighpathVersion,'.mat')), 'ff');
    ff = m.ff;

    % Load axon agglomerates
    m = load(fullfile(dataDir, sprintf('dendrites_%s.mat',suffixDendrites)));
    dendrites = m.dendrites(m.indBigDends);
    dendrites = arrayfun(@Agglo.fromSuperAgglo, dendrites, 'UniformOutput', false);
    clear m;

    % Get flight paths without comment but with start node
    ff = structfun(@(x)x(cellfun(@isempty, ff.comments)), ff, 'uni', 0);
    ff = structfun(@(x)x(~cellfun(@isempty, ff.startNode)), ff, 'uni', 0);

    segmentsLeftover = [];

    % Calculate overlap of all queries with segments
    [uniqueSegments, neighboursStartNode, ~, ~] = cellfun(@connectEM.queryAnalysis, ...
            ff.segIds, ff.neighbours, ff.nodes', ff.startNode', 'uni', 0);
    % Determine all overlaps of agglomerations with given queries
    [~, queryOverlap] = connectEM.queryAgglomerationOverlap(dendrites, segmentsLeftover, uniqueSegments, neighboursStartNode);
    % Make decision(s), here evidence/occurence threshold is applied
    % Always one (or none if evidence below 14, 1/2 node) start eqClass
    startAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 13), queryOverlap.start, 'uni', 0);
    % Exclude all queries that do not have a clear starting point
    idxNoClearStart = cellfun('isempty', startAgglo);
    % Multiple ends (all above 53vx evidence, corresponds to 2 full nodes)
    % Exclude startAgglo from endAgglo (as we do not want to count self-attachment)
    endAgglo = arrayfun(@(x) setdiff(queryOverlap.ends(x).eqClasses(queryOverlap.ends(x).occurences > 107),startAgglo{x}), (1:numel(queryOverlap.ends))', 'uni', 0);
    % Exclude all queries that do not have (at least one) clear end
    idxNoClearEnd = cellfun('isempty', endAgglo);
    % find multiple end agglo hits (possible if state of query agglos was
    % different from current agglo state
    idxMultiEdge = find(cellfun(@numel,endAgglo) > 1);
    % only those multi hits are problematic where interim agglos are not
    % fully covered by flight path (measured by checking if picked up
    % segments make up more than 85% of interim agglo volume)
    idxStrangeMultiEdge = idxMultiEdge( arrayfun(@(x)  sum(cellfun(@(y) sum(segmentMeta.voxelCount(y(ismember(y,ff.segIds{idxMultiEdge(x)})))),dendrites(endAgglo{idxMultiEdge(x)})) < cellfun(@(y) sum(segmentMeta.voxelCount(y)), dendrites(endAgglo{idxMultiEdge(x)}))*0.85) > 1, (1:numel(idxMultiEdge))')); % sum > 1 is because real end agglo is in most cases not full covered by flgiht path
    fprintf('%d multi hit flight paths detected, of which %d cannot be taken.\n',numel(idxMultiEdge),numel(idxStrangeMultiEdge))
    idxStrangeMultiEdge = accumarray(idxStrangeMultiEdge,1,[numel(idxNoClearEnd),1]);
    % 18.5% of queries excluded overall due to missing start or end (or both)
    idxGood = ~(idxNoClearStart | idxNoClearEnd | idxStrangeMultiEdge);
    % Display some statistics
    disp([num2str(sum(idxNoClearStart)./numel(idxNoClearStart)*100, '%.2f') '% of remaining queries have no clear start']);
    disp([num2str(sum(idxNoClearStart)) ' in total']);
    disp([num2str(sum(idxNoClearEnd)./numel(idxNoClearEnd)*100, '%.2f') '% of remaining queries have no clear end']);
    disp([num2str(sum(idxNoClearEnd)) ' in total']);
    disp([num2str(sum(idxGood)./numel(idxGood)*100, '%.2f') '% of remaining queries have clear start and ending']);
    disp([num2str(sum(idxGood)) ' in total']);
    disp([num2str(numel(cat(2, endAgglo{idxGood}))) ' attachments made by ' num2str(sum(idxGood)) ' queries']);
    % Find CC of eqClasses to be joined including single eqClasses with or
    % without dangling query
    edgesCC = cellfun(@(x,y)combnk([x y], 2), startAgglo(idxGood), endAgglo(idxGood), 'uni', 0);
    edgesCC = cat(1,edgesCC{:});
    edgesCC = sort(edgesCC, 2);
    eqClassCC = Graph.findConnectedComponents(edgesCC, true, true);
    sizeEqClassCC = sort(cellfun(@numel, eqClassCC), 'descend');
    eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(dendrites), cell2mat(eqClassCC)))'];
    display(sizeEqClassCC(1:10));

    results = struct;
    results.startAgglo = startAgglo;
    results.endAgglo = endAgglo;
    results.ff = ff;
    results.idxGood = idxGood;
    results.gitInfo = Util.gitInfo();

    % Save results and deprive writing permission
    saveFile = fullfile(dataDir, strcat('dendriteQueryOverlaps_',suffixVersion,'.mat'));
    save(saveFile, 'results', 'queryOverlap', 'idxNoClearStart', 'idxNoClearEnd');
    system(['chmod -w ' saveFile])
    
    saveFile = fullfile(dataDir, strcat('dendritePostQueryAnalysisState_',suffixVersion,'.mat'));
    save(saveFile);
    system(['chmod -w ' saveFile])

end

