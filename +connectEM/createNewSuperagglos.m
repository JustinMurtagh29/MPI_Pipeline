function createNewSuperagglos(param,state,casesToMerge)

    dataDir = fullfile(param.saveFolder, 'aggloState');

    [skeletonFolders, suffixFlightPaths, suffix, axonVersion] = connectEM.setQueryState(state);

    % Load current state of agglomerates
    superAgglos = load(fullfile(dataDir, strcat('axons_',num2str(axonVersion,'%.2i'),'.mat')));
    origAgglos = arrayfun(@Agglo.fromSuperAgglo, superAgglos.axons, 'uni', 0);
    superAgglos = superAgglos.axons(superAgglos.indBigAxons);

    load(fullfile(dataDir, strcat('caseDistinctions',suffix,'.mat')),'caseDistinctions',...
        'linkagesAgglos','flightPaths','linkagesFlat');
    load(fullfile(dataDir, strcat('attachedEndings',suffix,'.mat')),'executedFlightPaths');
    execusion = false(size(linkagesAgglos,1),1);
    execusion(executedFlightPaths) = true;
    executedFlightPaths = execusion;

    % Choose cases we use to merge agglos to super agglos
    if nargin < 3
        casesToMerge = [3 5 9 10];
    end
%     edgesCC = linkagesAgglos(ismember(caseDistinctions, casesToMerge),:);
    edgesCC = linkagesAgglos(executedFlightPaths,:);
    nans = isnan(edgesCC(:,1)) | isnan(edgesCC(:,2));
    edgesCC = edgesCC(~nans,:);
    edgesCC = sort(edgesCC, 2);
    eqClassCC = Graph.findConnectedComponents(edgesCC, true, true);
    eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(superAgglos), cell2mat(eqClassCC)))'];

    % Choose cases that are taken into account for super agglos
%     flightPaths.startAgglo(~ismember(caseDistinctions, casesToMerge),:) = [];
    flightPaths.startAgglo(~executedFlightPaths,:) = [];
%     flightPaths.endAgglo(~ismember(caseDistinctions, casesToMerge),:) = [];
    flightPaths.endAgglo(~executedFlightPaths,:) = [];
%     flightPaths.ff = structfun(@(x)x(ismember(caseDistinctions, casesToMerge)), flightPaths.ff, 'uni', 0);
    flightPaths.ff = structfun(@(x)x(executedFlightPaths), flightPaths.ff, 'uni', 0);
%     linkagesFlat = linkagesFlat(ismember(caseDistinctions, casesToMerge));
    linkagesFlat = linkagesFlat(executedFlightPaths);

    % Eliminate duplicate flight paths
    [~, positionUniques] = unique(sort(linkagesFlat,2), 'rows');
    duplicates = true(size(linkagesFlat,1),1);
    duplicates(positionUniques) = 0;

    flightPaths.startAgglo(duplicates,:) = [];
    flightPaths.endAgglo(duplicates,:) = [];
    flightPaths.ff = structfun(@(x)x(~duplicates), flightPaths.ff, 'uni', 0);

    axonsNew = connectEM.mergeSuperagglosBasedOnFlightPath(superAgglos, eqClassCCfull,...
        flightPaths.startAgglo, flightPaths.endAgglo, flightPaths.ff);
    axons = cat(1,axonsNew, superAgglos.axons(~superAgglos.indBigAxons));
    indBigAxons = false(length(axons));
    indBigAxons(1:length(axonsNew)) = true;

    % Save super agglos and deprive writing permission
    saveFile = fullfile(dataDir, strcat('axons_',num2str(axonVersion+1,'%.2i'),'.mat'));
    save(saveFile, 'axons','indBigAxons');
    system(['chmod -w ' saveFile])

end
