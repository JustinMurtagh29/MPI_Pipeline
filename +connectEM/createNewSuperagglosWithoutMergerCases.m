function createNewSuperagglosWithoutMergerCases(param,state,casesToMerge)

    dataDir = fullfile(param.saveFolder, 'aggloState');

    [skeletonFolders, suffixFlightPaths, suffix, axonVersion] = connectEM.setQueryState(state);

    % Load current state of agglomerates
    agglos = load(fullfile(dataDir, strcat('axons_',num2str(axonVersion,'%.2i'),'.mat')));
    origAgglos = arrayfun(@Agglo.fromSuperAgglo, agglos.axons, 'uni', 0);
    superAgglos = agglos.axons(agglos.indBigAxons);

    % Load linkages and cases for execusion
    load(fullfile(dataDir, strcat('caseDistinctions',suffix,'.mat')),...
        'linkagesAgglos','flightPaths','linkagesFlat');
    load(fullfile(dataDir, 'flightPathsWithoutMergerCases.mat'),'executedFlightPaths');
    execusion = false(size(linkagesAgglos,1),1);
    execusion(executedFlightPaths) = true;
    executedFlightPaths = execusion;

    % Find connected components of equivalent classes
    edgesCC = linkagesAgglos(executedFlightPaths,:);
    nans = isnan(edgesCC(:,1)) | isnan(edgesCC(:,2));
    edgesCC = edgesCC(~nans,:);
    edgesCC = sort(edgesCC, 2);
    eqClassCC = Graph.findConnectedComponents(edgesCC, true, true);
    % Generate new merged classes
    eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(superAgglos), cell2mat(eqClassCC)))'];

    % Use only those flight paths for execusion with the chosen case IDs
    flightPaths.startAgglo(~executedFlightPaths,:) = [];
    flightPaths.endAgglo(~executedFlightPaths,:) = [];
    flightPaths.ff = structfun(@(x)x(executedFlightPaths), flightPaths.ff, 'uni', 0);
    linkagesFlat = linkagesFlat(executedFlightPaths);

    % Eliminate duplicate flight paths
    [~, positionUniques] = unique(sort(linkagesFlat,2), 'rows');
    duplicates = true(size(linkagesFlat,1),1);
    duplicates(positionUniques) = 0;
    flightPaths.startAgglo(duplicates,:) = [];
    flightPaths.endAgglo(duplicates,:) = [];
    flightPaths.ff = structfun(@(x)x(~duplicates), flightPaths.ff, 'uni', 0);

    % Merging to super agglos
    axonsNew = connectEM.mergeSuperagglosBasedOnFlightPath(superAgglos, eqClassCCfull,...
        flightPaths.startAgglo, flightPaths.endAgglo, flightPaths.ff);
    % Concatenate small axons below 5 um
    axons = cat(1,axonsNew, agglos.axons(~agglos.indBigAxons));
    indBigAxons = false(length(axons),1);
    indBigAxons(1:length(axonsNew),1) = true;

    % Save super agglos and deprive writing permission
    saveFile = fullfile(dataDir, 'axons_05_withoutMergerCases.mat');
    save(saveFile, 'axons','indBigAxons');
    system(['chmod -w ' saveFile])

end
