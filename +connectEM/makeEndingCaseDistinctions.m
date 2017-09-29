function makeEndingCaseDistinctions(param,state)
    % makeEndingCaseDistinctions(param,state)
    %   Contrary to its name this function classifies flight paths and not
    %   endings. Flight paths which reach multiple agglomerates are ignored
    %   in this analysis.
    %
    % Written by
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    dataDir = fullfile(param.saveFolder, 'aggloState');
    
    [~, ~, suffix] = connectEM.setQueryState(state);

    endings = load(fullfile(dataDir, 'axonEndings.mat'));
    endingOverlap = load(fullfile(dataDir, strcat('axonEndingOverlaps',suffix,'.mat')));
    flightPaths = load(fullfile(dataDir, strcat('axonPostQueryAnalysisState',suffix,'.mat')), 'ff', 'startAgglo', 'endAgglo');
    
    % NOTE(amotta):
    % Q: Do multiple hits only occur for multiple endings from the same
    %    agglomerate or also for multiple endings from different
    %    agglomerates?
    % A: `connectEM.flightEndingOverlap` returns at most one ending per
    %    reached agglomerate. Multiple hits thus means that multiple
    %    agglomerates have overlap!
    
    % Which endings are linked together
    [linkages, idxMultipleHits] = cellfun(@connectEM.edgeCreator, endingOverlap.startEndingOverlaps, endingOverlap.endEndingOverlaps, 'uni', 0);
    linkagesFlat = cat(1, linkages{:});
    idxMultipleHits = cat(1, idxMultipleHits{:});
    
    % Exclude flight paths which reach **multiple agglomerates**.
    flightPaths.ff = structfun(@(x)x(~idxMultipleHits), flightPaths.ff, 'uni', 0);
    flightPaths.startAgglo(idxMultipleHits) = [];
    flightPaths.endAgglo(idxMultipleHits) = [];

    % Which agglos are linked together
    linkagesAgglos = cellfun(@connectEM.edgeCreator, flightPaths.startAgglo, flightPaths.endAgglo, 'uni', 0);
    linkagesAgglos = cat(1, linkagesAgglos{:});

    % Which endings belong to which agglo
    clusterSizes = cellfun(@max, endings.borderClusters);
    clusterLookup = repelem(1:numel(clusterSizes), clusterSizes);
    
    % NOTE(amotta): The ending indices in `linkages` and `linkagesFlat` can
    % take on the following valid values:
    % • nan ←→ no agglomerates was reached
    % • zero ←→ exactly one agglomerate, but no endings was reached
    % • positive integer ←→ exactly one ending was reached
    
    caseDistinctions = zeros(size(linkagesFlat,1),1);
    
    % Choose 1 for no attachment
    caseDistinctions(isnan(linkagesFlat(:,1)) & isnan(linkagesFlat(:,2))) = 1;
    
    % NOTE(amotta): Cases 2 and 3 do reach an agglomerate at exactly one
    % end of the flight path. However, we throw away the knowledge of
    % whether these paths also reached an ending or not!
    
    % Choose 2 for attachment at end but not at start
    caseDistinctions(isnan(linkagesFlat(:,1)) & ~isnan(linkagesFlat(:,2))) = 2;
    % Choose 3 for attachment to start but not at end
    caseDistinctions(~isnan(linkagesFlat(:,1)) & isnan(linkagesFlat(:,2))) = 3;
    
    % Find indices did match criteria so far
    % NOTE(amotta): Only look at flight paths which attach to
    % **agglomerates** at both the start and end of it.
    idxMatched = any(isnan(linkagesFlat),2);
    % Choose 4 for attachted to start not at ending and to end not at ending
    caseDistinctions(linkagesFlat(:,1) == 0 & linkagesFlat(:,2) == 0 & ~idxMatched) = 4;
    % Choose 5 for attachted to start at ending and to end not at ending
    % NOTE(amotta): Precise meaning of case 5 changes further down!
    caseDistinctions(linkagesFlat(:,1) ~= 0 & linkagesFlat(:,2) == 0 & ~idxMatched) = 5;
    % Choose 6 for attachted to start not at ending and to end at ending
    caseDistinctions(linkagesFlat(:,1) == 0 & linkagesFlat(:,2) ~= 0 & ~idxMatched) = 6;
    
    % Find indices did not match criteria so far
    % NOTE(amotta): Only look at flight paths which attach to **endings**
    % at both the start and end of it.
    idxMatched = idxMatched | any(linkagesFlat == 0,2);
    % Choose 7 for self attachment (within ending)
    caseDistinctions((linkagesFlat(:,1) == linkagesFlat(:,2)) & ~idxMatched) = 7;
    
    % Find indices did not match criteria so far
    % NOTE(amotta): Only look at flight paths which attach to **two
    % different endings** at the start and end of it.
    idxMatched = idxMatched | (linkagesFlat(:,1) == linkagesFlat(:,2));
    % Choose 8 for self attachment (within agglo)
    caseDistinctions((linkagesAgglos(:,1) == linkagesAgglos(:,2)) & ~idxMatched) = 8;
    
    % Find indices did not match criteria so far
    % NOTE(amotta): Only look at flight paths which attach to **two
    % different endings of two different agglomerates**.
    idxMatched = idxMatched | (linkagesAgglos(:,1) == linkagesAgglos(:,2));
    % Choose 9 for normal attachment
    % NOTE(amotta): Precise meaning of case 9 changes further down!
    caseDistinctions(~idxMatched) = 9;

    %% Define exclusion zone near dataset center
    options.border = [3000; -3000];
    borderNm = repmat(options.border, 1, 3);
    borderVoxel = round(bsxfun(@times, 1./param.raw.voxelSize, borderNm));
    bboxSmall = param.bbox + borderVoxel';

    % Check which flight path passes through exclusion zone of the dataset
    % NOTE(amotta): `true` for flight paths which never leave the core.
    idxWithin = cellfun(@(x)all(all(bsxfun(@gt, x, bboxSmall(:, 1)'),2) ...
        & all(bsxfun(@lt, x, bboxSmall(:, 2)'),2),1), flightPaths.ff.nodes);

    % Choose 10 for dangling edge that reaches EODS
    % NOTE(amotta): We seeded queries in the core of the dataset. Flight
    % paths which **leave the core and do not attach to an agglomerate**
    % are marked here. This case also captures flights which temporarily
    % leave the core but end within it.
    % NOTE(amotta): Precise meaning of case 10 changes further down!
    caseDistinctions(~idxWithin & caseDistinctions' == 3) = 10;
    
    %% Start to look for inconsistencies between sets of flight paths
    % NOTE(amotta): This analysis looks more at endings than at flight
    % paths. Shouldn't this be moved to `makeCaseDistinctionsOnEndings`?
    
    % NOTE(amotta): Both of the following checks are independent of the
    % flight paths directions. Consider the two case 9 `linkages` A → B and
    % A → C, where A, B, and C are ending indices. They will be re-
    % classified from case 9 to case 12. But the same happens for the
    % `linkages` A → B and C → A.
    
    % Check case 5 & 9 for inconsitencies and give them special case
    % Move all triangles in case 5 to case 11
    idx5 = caseDistinctions == 5;
    uniqueLinkages5 = unique(cat(2, linkagesFlat(idx5,1), linkagesAgglos(idx5,2)), 'rows');
    uniqueEndings5 = unique(uniqueLinkages5(:,1));
    count = histc(uniqueLinkages5(:,1), uniqueEndings5);
    
    % NOTE(amotta): `uniqueEndings5(count > 1)` is the list of all endings
    % which give rise to **multiple case 5 flight paths which reach
    % different agglomerates**.
    idxFlight5 = any(ismember(linkagesFlat, uniqueEndings5(count > 1)),2);
    caseDistinctions(idx5 & idxFlight5) = 11;

    % Move all triangles in case 9 to case 12
    idx9 = caseDistinctions == 9;
    uniqueLinkages9 = unique(sort(linkagesFlat(idx9,:),2), 'rows');
    uniqueEndings9 = unique(uniqueLinkages9(:));
    count = histc(uniqueLinkages9(:), uniqueEndings9);
    idxFlight9 = any(ismember(linkagesFlat, uniqueEndings9(count > 1)),2);
    caseDistinctions(idx9 & idxFlight9) = 12;
    
    %% Define exclusion zone for true ending definition
    options.border = [500; -500];
    borderNm = repmat(options.border, 1, 3);
    borderVoxel = round(bsxfun(@times, 1./param.raw.voxelSize, borderNm));
    bboxSmall = param.bbox + borderVoxel';
    
    % Check which flight path passes through exclusion zone of the dataset
    idxWithin = cellfun(@(x)all(all(bsxfun(@gt, x, bboxSmall(:, 1)'),2) ...
        & all(bsxfun(@lt, x, bboxSmall(:, 2)'),2),1), flightPaths.ff.nodes);

    % Choose 13 for definition of dangling edge with true ending
    caseDistinctions(~idxWithin & caseDistinctions' == 10) = 13;
    
    % NOTE(amotta): Case 10 now contains all flight paths which **cross the
    % 3 µm border, but the 500 nm one**. Case 13 contains the flight paths
    % which also reach the 500 nm border region. Interpretation unclear?

    %% Output count & frequency of each case
    tabulate(caseDistinctions)

    % Chose cases that are defined as attachted
    casesToCountAttached = [5 9 10];
    attachedQueries = linkagesFlat(ismember(caseDistinctions, casesToCountAttached),:);
    attachedEndings = unique(attachedQueries(:));

    % Save and deprive writing permission
    saveFile = fullfile(dataDir, strcat('caseDistinctions',suffix,'.mat'));
    save(saveFile);
    system(['chmod -w ' saveFile])

end

