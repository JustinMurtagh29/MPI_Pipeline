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
    endingCaseDistinctions(notProcessed) = 13;
    % Endings beyond 3um cutoff. Will be overwritten for those endings
    % attached by any constellation following below, so ouside without
    % attachment
    endingCaseDistinctions(outside) = 14;
    
    % Cases with redudancy, set those with more than 50% consensus as
    % solved
    for i=1:length(redundancies)
        if redundancies(i) == 1 & notProcessed(i) == 0
            tbl=tabulate(endingCases{i});
            if any(tbl(:,3) > 50)
                caseID = tbl(find(tbl(:,3) > 50),1);
                endingCaseDistinctionsMulti(i,1) = caseID;
            else
                endingCaseDistinctionsMulti(i,1) = NaN;
            end
        end
    end
    % Set redundancy cases
    % Start ending and end ending (with more than 50% consensus redundancy)
    endingCaseDistinctions(endingCaseDistinctionsMulti == 9) = 2;
    % Dangling flight path to EoDS 3um (with more than 50% consensus redundancy)
    endingCaseDistinctions(endingCaseDistinctionsMulti == 10) = 4;
    % Dangling flight path to EoDS 0.5um (with more than 50% consensus redundancy)
    endingCaseDistinctions(endingCaseDistinctionsMulti == 13) = 6;
    % Start ending and end agglo not an ending (with more than 50% consensus redundancy)
    endingCaseDistinctions(endingCaseDistinctionsMulti == 5) = 8;
    % Dangling flight path within DS (with more than 50% consensus redundancy)
    endingCaseDistinctions(endingCaseDistinctionsMulti == 3) = 10;
    % Contradicting redundancy
    endingCaseDistinctions(isnan(endingCaseDistinctionsMulti)) = 12;
    % Merger case E5a
    endingCaseDistinctions(endingCaseDistinctionsMulti == 12) = 15;
    % Merger case E5c
    endingCaseDistinctions(endingCaseDistinctionsMulti == 11) = 16;
    
    
    % Check the redundancy of dangling flight paths. For those which end
    % nodes have more than 1um distance another redundancy condition is
    % introduced, again min. 50% of the flight paths must have less than
    % the 1um distance.
    danglings = find(endingCaseDistinctions == 10);
    for i=1:length(danglings)
        danglingFlights = flightsOfEndingCases{danglings(i)}(endingCases{danglings(i)} == 3);
        endNodes = cellfun(@(x)x(end,:),ff.nodes(danglingFlights),'uni',0);
        distances = pdist(bsxfun(@times,cell2mat(endNodes'), [11.24 11.24 28]));
        if any(distances > 1000)
            if sum(distances < 1000) / numel(distances) <= 0.5
                endingCaseDistinctions(danglings(i)) = 11;
            end
        end
    end
        
    % E5 cases:
    casesE5 = find(endingCaseDistinctions == 12);
%     flightsE5 = flightsOfEndingCases(endingCaseDistinctions == 12);
    casesE5Distinctions = zeros(length(casesE5),1);
    endingCasesDanglingChecked = endingCases;
    flightsOfEndingCasesDanglingChecked = flightsOfEndingCases;
    for i=1:length(casesE5)
        % E5b
        danglings = flightsOfEndingCases{casesE5(i)}(find(endingCases{casesE5(i)} == 3));
        others = flightsOfEndingCases{casesE5(i)}(find(endingCases{casesE5(i)} ~= 3));
        endNodes = cellfun(@(x)x(end,:),ff.nodes(danglings),'uni',0);
        endNodes = cell2mat(endNodes');
        for j=1:length(danglings)
            distances = cellfun(@(x)min(pdist2(bsxfun(@times,endNodes, [11.24 11.24 28]),...
                bsxfun(@times,x, [11.24 11.24 28]))),ff.nodes(others));
            if any(distances < 500)
                endingCasesDanglingChecked{casesE5(i)}(j) = [];
                flightsOfEndingCasesDanglingChecked{casesE5(i)}(j) = [];
            end
            if numel(endingCasesDanglingChecked{casesE5(i)}) == 1
                endingCaseDistinctionsSingle(casesE5(i),1) = endingCasesDanglingChecked{casesE5(i)};
            end
        end
    end
    
    % Start ending and end ending (single flight path)
    endingCaseDistinctions(endingCaseDistinctionsSingle == 9) = 1;
    % Dangling flight path to EoDS 3um (single flight path)
    endingCaseDistinctions(endingCaseDistinctionsSingle == 10) = 3;
    % Dangling flight path to EoDS 0.5um (single flight path)
    endingCaseDistinctions(endingCaseDistinctionsSingle == 13) = 5;
    % Start ending and end agglo not an ending (single flight path)
    endingCaseDistinctions(endingCaseDistinctionsSingle == 5) = 7;
    % Dangling flight path within DS (single flight path)
    endingCaseDistinctions(endingCaseDistinctionsSingle == 3) = 9;
    % Merger case E5a
    endingCaseDistinctions(endingCaseDistinctionsSingle == 12) = 15;

        tabulate(endingCaseDistinctions)
        
%         endingCasesDanglingChecked
    casesE5 = find(endingCaseDistinctions == 12);
%     casesE5 = endingCases(endingCaseDistinctions == 12);
    casesE5Distinctions = zeros(length(casesE5),1);
       
    for i=1:length(casesE5)
        currentEnding = unique(endingCasesDanglingChecked{casesE5(i)});
        % E5b
        if isequal(currentEnding,[5 9])
            endingCaseDistinctions(casesE5(i),1) = 17;
        end
        % E5d
        if ismember(9,currentEnding)
            if any(ismember([3 10 13],currentEnding))
                endingCaseDistinctions(casesE5(i),1) = 18;
            end
        end
        % E5e
        if ismember(5,currentEnding)
            if any(ismember([3 10 13],currentEnding))
                endingCaseDistinctions(casesE5(i),1) = 19;
            end
        end
        if any(ismember([2 4 6],currentEnding))
            endingCaseDistinctions(casesE5(i),1) = 20;
        end
        % E5f
        if isequal(currentEnding,[3 10])
            endingCaseDistinctions(casesE5(i),1) = 21;
        end
        if isequal(currentEnding,[3 13])
            endingCaseDistinctions(casesE5(i),1) = 21;
        end
        if isequal(currentEnding,[10 13])
            endingCaseDistinctions(casesE5(i),1) = 21;
        end
        if isequal(currentEnding,[3 10 13])
            endingCaseDistinctions(casesE5(i),1) = 21;
        end
        if any(ismember([12],currentEnding))
            endingCaseDistinctions(casesE5(i),1) = 22;
        end
        
    end
     
    tabulate(endingCaseDistinctions)
    
    casesToMerge = [1:10];
    executedFlightPaths = flightsOfEndingCasesDanglingChecked(ismember(endingCaseDistinctions, casesToMerge));
    executedFlightPaths = cat(2,executedFlightPaths{:})';
   
    
    % Save and deprive writing permission
    saveFile = fullfile(dataDir, strcat('attachedEndings',suffix,'.mat'));
    save(saveFile,'attachedEndings');
    system(['chmod -w ' saveFile])
    
    tabulate(endingCaseDistinctions)
end
    