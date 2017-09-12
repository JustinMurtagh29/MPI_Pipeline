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
    outside = ~inside;
    
    endingOverlap = load(fullfile(dataDir, strcat('axonEndingOverlaps',suffix,'.mat')));
    load(fullfile(dataDir, strcat('caseDistinctions',suffix,'.mat')),'caseDistinctions','idxMultipleHits');
    
    startEndingOverlaps=endingOverlap.startEndingOverlaps;
    startEndingOverlaps = startEndingOverlaps(~idxMultipleHits);
    endEndingOverlaps=endingOverlap.endEndingOverlaps;
    endEndingOverlaps = endEndingOverlaps(~idxMultipleHits);
    
    clusterSizes = cellfun(@max, endingData.borderClusters);
    
    endingCases = cell(sum(clusterSizes),1);
    for i=1:length(caseDistinctions)
        if startEndingOverlaps{i} > 0
            endingCases{startEndingOverlaps{i},1} = [endingCases{startEndingOverlaps{i},1} caseDistinctions(i)];
        end
    end
    
    for i=1:length(caseDistinctions)
        if numel(endEndingOverlaps{i}) < 2
            if endEndingOverlaps{i} > 0
                endingCases{endEndingOverlaps{i},1} = [endingCases{endEndingOverlaps{i},1} caseDistinctions(i)];
            end
        end
    end
    
    endingCaseDistinctionsSingle = zeros(size(endingCases,1),1);
    notProcessed = cellfun('isempty',endingCases);
    endingCaseDistinctionsSingle(notProcessed) = 16;
    endingCaseDistinctionsSingle(outside) = 15;
    endingCaseDistinctionsMulti = endingCaseDistinctionsSingle;
    
    
    redundancy = cellfun(@numel,endingCases);
    redundancies = redundancy > 1;
    % cases without redanduncy
    for i=1:length(redundancies)
        if redundancies(i) == 0 & notProcessed(i) == 0
            endingCaseDistinctionsSingle(i,1) = endingCases{i};
        end
    end
    
    endingCaseDistinctions = zeros(size(endingCases,1),1);
    endingCaseDistinctions(notProcessed) = 12;
    endingCaseDistinctions(outside) = 13;
    
    endingCaseDistinctions(endingCaseDistinctionsSingle == 9) = 1;
    endingCaseDistinctions(endingCaseDistinctionsSingle == 10) = 3;
    endingCaseDistinctions(endingCaseDistinctionsSingle == 13) = 5;
    endingCaseDistinctions(endingCaseDistinctionsSingle == 5) = 7;
    endingCaseDistinctions(endingCaseDistinctionsSingle == 3) = 9;
    
    % cases with redudancy
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
    endingCaseDistinctions(endingCaseDistinctionsMulti == 9) = 2;
    endingCaseDistinctions(endingCaseDistinctionsMulti == 10) = 4;
    endingCaseDistinctions(endingCaseDistinctionsMulti == 13) = 6;
    endingCaseDistinctions(endingCaseDistinctionsMulti == 5) = 8;
    endingCaseDistinctions(endingCaseDistinctionsMulti == 3) = 10;
    endingCaseDistinctions(isnan(endingCaseDistinctionsMulti)) = 11;
    
    % Chose cases that are defined as attachted
    casesToCountAttached = [1 2 3 4 5 6 7 8 10];
    attachedEndings = find(ismember(endingCaseDistinctions, casesToCountAttached),:);
    
     % Save and deprive writing permission
    saveFile = fullfile(dataDir, strcat('attachedEndings',suffix,'.mat'));
    save(saveFile);
    system(['chmod -w ' saveFile])


    
    tabulate(endingCaseDistinctions)
end
    