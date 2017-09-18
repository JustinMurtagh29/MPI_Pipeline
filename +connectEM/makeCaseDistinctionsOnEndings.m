function makeCaseDistinctionsOnEndings(param,state)
    
    dataDir = fullfile(param.saveFolder, 'aggloState');
    
    [skeletonFolders, suffixFlightPaths, suffix] = connectEM.setQueryState(state);    

    % Load data from ending generation
    endingData = fullfile(dataDir, 'axonEndingsAllData.mat');
    endingData = load(endingData);
    % Extract some variables
    idxAll = endingData.idxAll;
    idxCanidateFound = endingData.axonMask;
    borderIds = endingData.borderIds;
    borderClusters = endingData.borderClusters ;

    % Load directionality information
    directionality = fullfile(dataDir, 'axonEndingInputData.mat');
    directionality = load(directionality, 'directionality');
    directionality = directionality.directionality;

    % Load border CoMs
    borderCoM = fullfile(param.saveFolder, 'globalBorder.mat');
    borderCoM = load(borderCoM, 'borderCoM');
    borderCoM = borderCoM.borderCoM;
    borderPositions = cellfun(@(x) borderCoM(x,:), directionality.borderIdx(idxCanidateFound),'uni',0);

    % Find indices of ending candidates in directionality lists (scores,
    % pca...) and exclude redundant endings
    idxUse = idxAll(idxCanidateFound);
    idxCluster = [];
    for j=1:length(borderClusters)
        idxCluster{j,1} = [];
        for k=1:max(borderClusters{j})
            idxCluster{j,1}{k,1} = idxUse{j,1}(find(borderClusters{j}==k));
        end
    end

    % Write out absolut values of seg direction scores
    SegDirScores = cellfun(@(x)abs(x),directionality.scores(idxCanidateFound),'uni',0);

    % Dataset border conditions
    options.border = [3000; -3000];
    borderNm = repmat(options.border, 1, 3);
    borderVoxel = round(bsxfun(@times, 1./param.raw.voxelSize, borderNm));
    bboxSmall = param.bbox + borderVoxel';

    % Keep only candidate with highest score for each cluster and exclude
    % candidate outside border cutoff
    candidateUse = {};
    for j=1:length(idxCluster)
        candidateUse{j,1} = [];
        for k=1:length(idxCluster{j})
            scoreClusters = SegDirScores{j,1}(idxCluster{j,1}{k,1});
            [~, sortIdx] = sort(scoreClusters, 'descend');
            candidate = idxCluster{j,1}{k,1}(sortIdx(1));
            if all(bsxfun(@gt, borderPositions{j}(candidate,:), bboxSmall(:, 1)'),2)...
                    & all(bsxfun(@lt, borderPositions{j}(candidate,:), bboxSmall(:, 2)'),2)
                candidateUse{j,1}(end+1,1) = candidate;
            else
                candidateUse{j,1}(end+1,1) = 0;
            end
        end
    end

    
    inside = cat(1,candidateUse{:});
    inside = inside > 0;
    % Logical matrix which indicates endings beyond cutoff border
    outside = ~inside;
    
    % Load flight path - ending overlap, case distinctions on flight
    % paths and flight paths itself
    endingOverlap = load(fullfile(dataDir, strcat('axonEndingOverlaps',suffix,'.mat')));
    load(fullfile(dataDir, strcat('caseDistinctions',suffix,'.mat')),'caseDistinctions','idxMultipleHits','flightPaths');
    ff = flightPaths.ff;
    
    startEndingOverlaps = endingOverlap.startEndingOverlaps;
    startEndingOverlaps = startEndingOverlaps(~idxMultipleHits);
    endEndingOverlaps = endingOverlap.endEndingOverlaps;
    endEndingOverlaps = endEndingOverlaps(~idxMultipleHits);
    
    clusterSizes = cellfun(@max, endingData.borderClusters);
    
    %% Case distinctions:
    % Read out caseID for all attachments at start ending 
    endingCases = cell(sum(clusterSizes),1);
    for i=1:length(caseDistinctions)
        if startEndingOverlaps{i} > 0
            endingCases{startEndingOverlaps{i},1} = [endingCases{startEndingOverlaps{i},1} caseDistinctions(i)];
        end
    end
    
    % Read out caseID for all attachments at end ending
    for i=1:length(caseDistinctions)
        if endEndingOverlaps{i} > 0
            endingCases{endEndingOverlaps{i},1} = [endingCases{endEndingOverlaps{i},1} caseDistinctions(i)];
        end
    end
    
    % Keep flight path ID for later use
    flightsOfEndingCases = cell(sum(clusterSizes),1);
    for i=1:length(caseDistinctions)
        if startEndingOverlaps{i} > 0
            flightsOfEndingCases{startEndingOverlaps{i},1} = [flightsOfEndingCases{startEndingOverlaps{i},1} i];
        end
    end
    for i=1:length(caseDistinctions)
        if endEndingOverlaps{i} > 0
            flightsOfEndingCases{endEndingOverlaps{i},1} = [flightsOfEndingCases{endEndingOverlaps{i},1} i];
        end
    end
    
    
    endingCaseDistinctionsSingle = zeros(size(endingCases,1),1);
    notProcessed = cellfun('isempty',endingCases);
    endingCaseDistinctionsMulti = endingCaseDistinctionsSingle;
    
    % Check redundancy
    redundancy = cellfun(@numel,endingCases);
    redundancies = redundancy > 1;
    % Cases without redundancy
    for i=1:length(redundancies)
        if redundancies(i) == 0 & notProcessed(i) == 0
            endingCaseDistinctionsSingle(i,1) = endingCases{i};
        end
    end
    
    % Set cases first for those without redundancy
    endingCaseDistinctions = zeros(size(endingCases,1),1);
    % Endings inside which are not attached by any flight path (beyond 3um
    % cutoff will be overwritten in next step)
    endingCaseDistinctions(notProcessed) = 15;
    % Endings beyond 3um cutoff. Will be overwritten for those endings
    % attached by any constellation following below, so ouside without
    % attachment
    endingCaseDistinctions(outside) = 16;
    
    % Eliminiate cases without start attachment at endings (only for
    % redundant endings since the old queries disturb the analysis at this
    % point.
    redundantCases = find(redundancies == 1);
    for i=1:length(redundantCases)
        noStartAttachments = cell2mat(arrayfun(@(x)ismember(x,[2 4 6]),endingCases{redundantCases(i)},'uni',0));
        endingCases{redundantCases(i)}(noStartAttachments) = [];
        flightsOfEndingCases{redundantCases(i)}(noStartAttachments) = [];
        if numel(endingCases{redundantCases(i)}) == 1
            endingCaseDistinctionsSingle(redundantCases(i),1) = endingCases{redundantCases(i)};
        end
    end
    
    % Check redundancy of dangling flight paths, those which are getting 
    % continued by another flight path are no longer taken into account.   
    redundancy = cellfun(@numel,endingCases);
    redundancies = redundancy > 1;
    redundantCases = find(redundancies == 1);
    endingCasesDanglingChecked = endingCases;
    flightsOfEndingCasesDanglingChecked = flightsOfEndingCases;
    for i=1:length(redundantCases)
        danglings = flightsOfEndingCases{redundantCases(i)}(find(endingCases{redundantCases(i)} == 3));
        others = flightsOfEndingCases{redundantCases(i)}(find(endingCases{redundantCases(i)} ~= 3));
        endNodes = cellfun(@(x)x(end,:),ff.nodes(danglings),'uni',0);
        endNodes = cell2mat(endNodes');
        distances = [];
        if ~isempty(others)
            for j=1:length(danglings)
                distance = cell2mat(cellfun(@(x)min(pdist2(bsxfun(@times,endNodes(j,:), [11.24 11.24 28]),...
                    bsxfun(@times,x, [11.24 11.24 28]))),ff.nodes(others),'uni',0));
                distances(j,1) = min(distance);
            end
        end
        redundantDanglings = distances < 500;
        endingCasesDanglingChecked{redundantCases(i)}(redundantDanglings) = [];
        flightsOfEndingCasesDanglingChecked{redundantCases(i)}(redundantDanglings) = [];
        if numel(endingCasesDanglingChecked{redundantCases(i)}) == 1
            endingCaseDistinctionsSingle(redundantCases(i),1) = endingCasesDanglingChecked{redundantCases(i)};
        end
    end
    endingCases = endingCasesDanglingChecked;
    flightsOfEndingCases = flightsOfEndingCasesDanglingChecked;
    
    
    %% Set obvious cases for endings with just a single flight path:
    % Start ending and end ending (single flight path)
    endingCaseDistinctions(endingCaseDistinctionsSingle == 9) = 1;
    % Dangling flight path to EoDS 3um (single flight path)
    endingCaseDistinctions(endingCaseDistinctionsSingle == 10) = 3;
    % Dangling flight path to EoDS 0.5um (single flight path)
    endingCaseDistinctions(endingCaseDistinctionsSingle == 13) = 3;
    % Start ending and end agglo not an ending (single flight path)
    endingCaseDistinctions(endingCaseDistinctionsSingle == 5) = 5;
    % Dangling flight path within DS (single flight path)
    endingCaseDistinctions(endingCaseDistinctionsSingle == 3) = 7;
    % Merger case E5a
    endingCaseDistinctions(endingCaseDistinctionsSingle == 12) = 9;

    %% Cases with redudancy, set those with more than 50% consensus as
    % solved
    redundancy = cellfun(@numel,endingCases);
    redundancies = redundancy > 1;
    redundantCases = find(redundancies == 1);
    for i=1:length(redundantCases)
        tbl=tabulate(endingCases{redundantCases(i)});
        if any(tbl(:,3) > 50)
            caseID = tbl(find(tbl(:,3) > 50),1);
            endingCaseDistinctionsMulti(redundantCases(i),1) = caseID;
            otherCases = unique(endingCases{redundantCases(i)}(endingCases{redundantCases(i)} ~= caseID));
            if ~isempty(otherCases)
                for j=1:length(otherCases)
                    endingCases{redundantCases(i)}(endingCases{redundantCases(i)} == otherCases(j)) = [];
                    flightsOfEndingCases{redundantCases(i)}(endingCases{redundantCases(i)} == otherCases(j)) = [];
                end
            end
        else
            endingCaseDistinctionsMulti(redundantCases(i),1) = NaN;
        end
    end
    % Set redundancy cases
    % Start ending and end ending (with more than 50% consensus redundancy)
    endingCaseDistinctions(endingCaseDistinctionsMulti == 9) = 2;
    % Dangling flight path to EoDS 3um (with more than 50% consensus redundancy)
    endingCaseDistinctions(endingCaseDistinctionsMulti == 10) = 4;
    % Dangling flight path to EoDS 0.5um (with more than 50% consensus redundancy)
    endingCaseDistinctions(endingCaseDistinctionsMulti == 13) = 4;
    % Start ending and end agglo not an ending (with more than 50% consensus redundancy)
    endingCaseDistinctions(endingCaseDistinctionsMulti == 5) = 6;
    % Dangling flight path within DS (with more than 50% consensus redundancy)
    endingCaseDistinctions(endingCaseDistinctionsMulti == 3) = 8;
    % Contradicting redundancy
    endingCaseDistinctions(isnan(endingCaseDistinctionsMulti)) = 17;
    % Merger case E5a
    endingCaseDistinctions(endingCaseDistinctionsMulti == 12) = 9;
    % Merger case E5c
    endingCaseDistinctions(endingCaseDistinctionsMulti == 11) = 11;
    
    
    % Check the redundancy of dangling flight paths. For those which end
    % nodes have more than 1um distance another redundancy condition is
    % introduced, again min. 50% of the flight paths must have less than
    % the 1um distance.
    danglings = find(endingCaseDistinctions == 8);
    for i=1:length(danglings)
        danglingFlights = flightsOfEndingCases{danglings(i)}(endingCases{danglings(i)} == 3);
        endNodes = cellfun(@(x)x(end,:),ff.nodes(danglingFlights),'uni',0);
        distances = pdist(bsxfun(@times,cell2mat(endNodes'), [11.24 11.24 28]));
        if any(distances > 1000)
            if sum(distances < 1000) / numel(distances) <= 0.5
                endingCaseDistinctions(danglings(i)) = 18;
            end
        end
    end
        
    
    
        tabulate(endingCaseDistinctions)
        
    %% Remaining E5 cases:

    casesE5 = find(endingCaseDistinctions == 17);
       
    for i=1:length(casesE5)
        currentEnding = unique(endingCases{casesE5(i)});
        % E5b
        if isequal(currentEnding,[5 9])
            endingCaseDistinctions(casesE5(i),1) = 10;
        end
        % E5d
        if ismember(9,currentEnding)
            if any(ismember([3 10 13],currentEnding))
                endingCaseDistinctions(casesE5(i),1) = 12;
            end
        end
        % E5e
        if ismember(5,currentEnding)
            if any(ismember([3 10 13],currentEnding))
                endingCaseDistinctions(casesE5(i),1) = 13;
            end
        end
        % E5f
        if ismember(3,currentEnding)
            if any(ismember([10 13],currentEnding))
                endingCaseDistinctions(casesE5(i),1) = 14;
            end
        end
        
        % Rest
        if isequal(currentEnding,[3 10])
            endingCaseDistinctions(casesE5(i),1) = 18;
        end
        if isequal(currentEnding,[3 13])
            endingCaseDistinctions(casesE5(i),1) = 18;
        end
        if isequal(currentEnding,[3 10 13])
            endingCaseDistinctions(casesE5(i),1) = 18;
        end
        
        if isequal(currentEnding,[10 13])
            endingCaseDistinctions(casesE5(i),1) = 18;
        end
        
        if any(ismember([12],currentEnding))
            endingCaseDistinctions(casesE5(i),1) = 19;
        end
        
    end
    
    % Check distance between all remaining dangling flight paths. With more
    % than 2 um add to case E5 f
    danglings = find(endingCaseDistinctions == 18);
    for i=1:length(danglings)
        danglingFlights = flightsOfEndingCases{danglings(i)}...
            (endingCases{danglings(i)} == 3 | endingCases{danglings(i)} == 10 | endingCases{danglings(i)} == 13);
        endNodes = cellfun(@(x)x(end,:),ff.nodes(danglingFlights),'uni',0);
        distances = pdist(bsxfun(@times,cell2mat(endNodes'), [11.24 11.24 28]));
        if any(distances > 2000)
            endingCaseDistinctions(danglings(i)) = 14;
        end
    end
    
    % Sum up all remaining configurations
    endingCaseDistinctions(endingCaseDistinctions == 0) = 17;
    endingCaseDistinctions(endingCaseDistinctions == 18) = 17;
    endingCaseDistinctions(endingCaseDistinctions == 19) = 17;
    
     
    tabulate(endingCaseDistinctions)
    
    % Choose cases for merging and generation of queries
    casesToMerge = [1:6, 8:14];
    executedFlightPaths = flightsOfEndingCasesDanglingChecked(ismember(endingCaseDistinctions, casesToMerge));
    executedFlightPaths = unique(cat(2,executedFlightPaths{:})');
   
    casesToCountAttached = [1:6 8:14 16 17];
    attachedEndings = find(ismember(endingCaseDistinctions, casesToCountAttached));

    
    
    % Save and deprive writing permission
    saveFile = fullfile(dataDir, strcat('attachedEndings',suffix,'.mat'));
    save(saveFile,'attachedEndings','executedFlightPaths');
    system(['chmod -w ' saveFile])
    
    tabulate(endingCaseDistinctions)
end
    