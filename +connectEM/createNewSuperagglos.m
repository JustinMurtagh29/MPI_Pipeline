function createNewSuperagglos(param,state,casesToMerge)

    dataDir = fullfile(param.saveFolder, 'aggloState');

    [skeletonFolders, suffixFlightPaths, suffix, axonVersion] = connectEM.setQueryState(state);

    % Load current state of agglomerates
    superAgglos = load(fullfile(dataDir, strcat('axons_',num2str(axonVersion,'%.2i'),'.mat')));
    origAgglos = arrayfun(@Agglo.fromSuperAgglo, superAgglos.axons, 'uni', 0);
    superAgglos = superAgglos.axons(superAgglos.indBigAxons);

    load(fullfile(dataDir, strcat('caseDistinctions',suffix,'.mat')),'caseDistinctions',...
        'linkagesAgglos','flightPaths','linkagesFlat');

    % Choose cases we use to merge agglos to super agglos
    if nargin < 3
        casesToMerge = [3 5 9 10];
    end
    edgesCC = linkagesAgglos(ismember(caseDistinctions, casesToMerge),:);
    nans = isnan(edgesCC(:,1)) | isnan(edgesCC(:,2));
    edgesCC = edgesCC(~nans,:);
    edgesCC = sort(edgesCC, 2);
    eqClassCC = Graph.findConnectedComponents(edgesCC, true, true);
    eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(superAgglos), cell2mat(eqClassCC)))'];

    % Choose cases that are taken into account for super agglos
    flightPaths.startAgglo(~ismember(caseDistinctions, casesToMerge),:) = [];
    flightPaths.endAgglo(~ismember(caseDistinctions, casesToMerge),:) = [];
    flightPaths.ff = structfun(@(x)x(ismember(caseDistinctions, casesToMerge)), flightPaths.ff, 'uni', 0);
    linkagesFlat = linkagesFlat(ismember(caseDistinctions, casesToMerge));

    % Eliminate duplicate flight paths
    [~, positionUniques] = unique(sort(linkagesFlat,2), 'rows');
    duplicates = true(size(linkagesFlat,1),1);
    duplicates(positionUniques) = 0;

    flightPaths.startAgglo(duplicates,:) = [];
    flightPaths.endAgglo(duplicates,:) = [];
    flightPaths.ff = structfun(@(x)x(~duplicates), flightPaths.ff, 'uni', 0);

    axons = connectEM.mergeSuperagglosBasedOnFlightPath(superAgglos, eqClassCCfull,...
        flightPaths.startAgglo, flightPaths.endAgglo, flightPaths.ff);

    % Save super agglos and deprive writing permission
    saveFile = fullfile(dataDir, strcat('axons_',num2str(axonVersion+1,'%.2i'),'.mat'));
    save(saveFile, 'axons');
    system(['chmod -w ' saveFile])

end
