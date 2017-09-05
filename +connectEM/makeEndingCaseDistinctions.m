function makeEndingCaseDistinctions(param)

    dataDir = fullfile(param.saveFolder, 'aggloState');

    superAgglos = load(fullfile(dataDir, 'axons_04.mat'));
    origAgglos = arrayfun(@Agglo.fromSuperAgglo, superAgglos.axons, 'uni', 0);
    superAgglos = superAgglos.axons(superAgglos.indBigAxons);

    endings = load(fullfile(dataDir, 'axonEndings.mat'));
    endingOverlap = load(fullfile(dataDir, 'axonEndingOverlaps.mat'));
    flightPaths = load(fullfile(dataDir, 'axonPostQueryAnalysisState.mat'), 'ff', 'startAgglo', 'endAgglo');
    clusterSizes = cellfun(@max, endings.borderClusters);

    % Which endings are linked together
    [linkages, idxMultipleHits] = cellfun(@(x,y)connectEM.edgeCreator(x,y), endingOverlap.startEndingOverlaps, endingOverlap.endEndingOverlaps, 'uni', 0);
    linkagesFlat = cat(1, linkages{:});
    idxMultipleHits = cat(1, idxMultipleHits{:});
    flightPaths.ff = structfun(@(x)x(~idxMultipleHits), flightPaths.ff, 'uni', 0);
    flightPaths.startAgglo(idxMultipleHits) = [];
    flightPaths.endAgglo(idxMultipleHits) = [];

    % Which agglos are linked together
    linkagesAgglos = cellfun(@(x,y)connectEM.edgeCreator(x,y), flightPaths.startAgglo, flightPaths.endAgglo, 'uni', 0);
    linkagesAgglos = cat(1, linkagesAgglos{:});

    % Which endings belong to which agglo
    clusterSizes = cellfun(@max, endings.borderClusters);
    clusterLookup = repelem(1:numel(clusterSizes), clusterSizes);

    caseDistinctions = zeros(size(linkagesFlat,1),1);
    % Choose 1 for no attachment
    caseDistinctions(isnan(linkagesFlat(:,1)) & isnan(linkagesFlat(:,2))) = 1;
    % Choose 2 for attachment at end but not at start
    caseDistinctions(isnan(linkagesFlat(:,1)) & ~isnan(linkagesFlat(:,2))) = 2;
    % Choose 3 for attachment to start but not at end
    caseDistinctions(~isnan(linkagesFlat(:,1)) & isnan(linkagesFlat(:,2))) = 3;
    % Find indices did match criteria so far
    idxMatched = any(isnan(linkagesFlat),2);
    % Choose 4 for attachted to start not at ending and to end not at ending
    caseDistinctions(linkagesFlat(:,1) == 0 & linkagesFlat(:,2) == 0 & ~idxMatched) = 4;
    % Choose 5 for attachted to start at ending and to end not at ending
    caseDistinctions(linkagesFlat(:,1) ~= 0 & linkagesFlat(:,2) == 0 & ~idxMatched) = 5;
    % Choose 6 for attachted to start not at ending and to end at ending
    caseDistinctions(linkagesFlat(:,1) == 0 & linkagesFlat(:,2) ~= 0 & ~idxMatched) = 6;
    % Find indices did not match criteria so far
    idxMatched = idxMatched | any(linkagesFlat == 0,2);
    % Choose 7 for self attachment (within ending)
    caseDistinctions((linkagesFlat(:,1) == linkagesFlat(:,2)) & ~idxMatched) = 7;
    % Find indices did not match criteria so far
    idxMatched = idxMatched | (linkagesFlat(:,1) == linkagesFlat(:,2));
    % Choose 8 for self attachment (within agglo)
    caseDistinctions((linkagesAgglos(:,1) == linkagesAgglos(:,2)) & ~idxMatched) = 8;
    % Find indices did not match criteria so far
    idxMatched = idxMatched | (linkagesAgglos(:,1) == linkagesAgglos(:,2));
    % Choose 9 for normal attachment
    caseDistinctions(~idxMatched) = 9;

    % Define exclusion zone near dataset center
    options.border = [3000; -3000];
    borderNm = repmat(options.border, 1, 3);
    borderVoxel = round(bsxfun(@times, 1./param.raw.voxelSize, borderNm));
    bboxSmall = param.bbox + borderVoxel';

    % Check which flight path passes through exclusion zone of the dataset
    idxWithin = cellfun(@(x)all(all(bsxfun(@gt, x, bboxSmall(:, 1)'),2) ...
        & all(bsxfun(@lt, x, bboxSmall(:, 2)'),2),1), flightPaths.ff.nodes);

    % Choose 10 for dangling edge that reaches EODS
    caseDistinctions(~idxWithin & caseDistinctions' == 3) = 10;

    % Check case 5 & 9 for inconsitencies and give them special case
    % Move all triangles in case 5 to case 11
    idx5 = caseDistinctions == 5;
    uniqueLinkages5 = unique(cat(2, linkagesFlat(idx5,1), linkagesAgglos(idx5,2)), 'rows');
    uniqueEndings5 = unique(uniqueLinkages5(:,1));
    count = histc(uniqueLinkages5(:,1), uniqueEndings5);
    idxFlight5 = any(ismember(linkagesFlat, uniqueEndings5(count > 1)),2);
    caseDistinctions(idx5 & idxFlight5) = 11;

    % Move all triangles in case 9 to case 12
    idx9 = caseDistinctions == 9;
    uniqueLinkages9 = unique(sort(linkagesFlat(caseDistinctions == 9,:),2), 'rows');
    uniqueEndings9 = unique(uniqueLinkages9(:));
    count = histc(uniqueLinkages9(:), uniqueEndings9);
    idxFlight9 = any(ismember(linkagesFlat, uniqueEndings9(count > 1)),2);
    caseDistinctions(idx9 & idxFlight9) = 12;

    % Output count & frequency of each case
    tabulate(caseDistinctions)

    % Chose cases that are defined as attachted
    casesToCountAttached = [5 9 10];
    attachedQueries = linkagesFlat(ismember(caseDistinctions, casesToCountAttached),:);
    attachedEndings = unique(attachedQueries(:));

    % Save
    save(fullfile(dataDir, 'attachedEndings.mat'));

    % Choose cases we use to merge agglos to super agglos
    casesToCountAttached = [5 9];
    edgesCC = linkagesAgglos(ismember(caseDistinctions, casesToCountAttached),:);
    edgesCC = sort(edgesCC, 2);
    eqClassCC = Graph.findConnectedComponents(edgesCC, true, true);
    eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(superAgglos), cell2mat(eqClassCC)))'];
    
    % Choose cases that are additionally taken into account for super agglos
    casesToCountAttached = [3 5 9 10];
    startAgglo(~ismember(caseDistinctions, casesToCountAttached),:) = [];
    endAgglo(~ismember(caseDistinctions, casesToCountAttached),:) = [];

    [uniqueLinkages, positionUniques] = unique(sort(linkagesFlat,2), 'rows');
    dublicates = ones(size(linkagesFlat,1),1);
    dublicates(positionUniques) = 0;
    dublicates = logical(dublicates);
%     uniqueLinkages2 = linkagesFlat(~dublicates,:);
%     uniqueLinkages3 = unique(sort(uniqueLinkages2,2), 'rows');
    
    startAgglo(dublicates,:) = [];
    endAgglo(dublicates,:) = [];
    
    axons = mergeSuperagglosBasedOnFlightPath(superAgglos, eqClassCCfull, startAgglo, endAgglo, ff)

    % Save super agglos
    save(fullfile(dataDir, 'axons_05.mat'));

end

