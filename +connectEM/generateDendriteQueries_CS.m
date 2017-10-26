function generateDendriteQueries_CS(param,suffix,queryRun)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>

    if ~exist('suffix','var')
        suffix = '';
    end

    dataDir = fullfile(param.saveFolder, 'aggloState');

    % Load data from ending generation
    endingData = fullfile(dataDir, sprintf('dendriteEndingsAllData%s.mat',suffix));
    endingData = load(endingData);
    % Extract some variables
    idxAll = endingData.idxAll;
    idxCanidateFound = endingData.dendriteMask;
    borderIds = endingData.borderIds;
    borderClusters = endingData.borderClusters ;

    % Load directionality information
    directionality = fullfile(dataDir, sprintf('dendriteEndingInputData%s.mat',suffix));
    directionality = load(directionality, 'directionality');
    directionality = directionality.directionality;

    % Load border CoMs
    borderCoM = fullfile(param.saveFolder, 'globalBorder.mat');
    borderCoM = load(borderCoM, 'borderCoM');
    borderCoM = borderCoM.borderCoM;
    borderPositions = cellfun(@(x) borderCoM(x,:), directionality.borderIdx(idxCanidateFound),'uni',0);

    % Load larger 5 micron agglomerates
    m = load(fullfile(dataDir, 'dendrites_03_v2_splitmerged.mat'));
    dendrites = m.dendrites(m.indBigDends);
    dendrites = arrayfun(@Agglo.fromSuperAgglo, dendrites, 'UniformOutput', false);
    display([num2str(numel(m.indBigDends)) ' dendrites in total']);
    display([num2str(sum(m.indBigDends)) ' larger 5um']);

    % Determine endings which are not redundant(already attached by flight path)
%     m = load(fullfile(dataDir, strcat('attachedEndings',suffix,'.mat')), 'endingCaseDistinctions');
%     casesToCountAttached = [1:4 6:17];
%     attachedEndings = find(ismember(m.endingCaseDistinctions, casesToCountAttached));
%     clear m
    attachedEndings = [];

    % Find indices of ending candidates in directionality lists (scores,
    % pca...) and exclude redundant endings
    idxUse = idxAll(idxCanidateFound);
    counter = 1;
    idxCluster = [];
    for j=1:length(borderClusters)
        idxCluster{j,1} = {};
        for k=1:max(borderClusters{j})
            if ismember(counter,attachedEndings)
                idxCluster{j,1}{k,1} = [];
            else
                idxCluster{j,1}{k,1} = idxUse{j,1}(find(borderClusters{j}==k));
            end
            counter = counter + 1;
        end
        idxCluster{j,1} = idxCluster{j,1}(~cellfun('isempty',idxCluster{j,1}));
    end

    % Exclude axons without a single ending left
    redundant = cellfun('isempty',idxCluster);
    idxCluster = idxCluster(~redundant);

    % Write out absolut values of seg direction scores
    SegDirScores = cellfun(@(x)abs(x),directionality.scores(idxCanidateFound),'uni',0);
    SegDirScores = SegDirScores(~redundant);
    borderPositions = borderPositions(~redundant);

    % Dataset border conditions
    options.border = [3000; -3000];
    borderNm = repmat(options.border, 1, 3);
    borderVoxel = round(bsxfun(@times, 1./param.raw.voxelSize, borderNm));
    bboxSmall = param.bbox + borderVoxel';

    % Keep only candidate with highest score for each cluster and exclude
    % candidate outside border cutoff
    candidateUse = {};
    candidateOutside = {};
    for j=1:length(idxCluster)
        candidateUse{j,1} = [];
        candidateOutside{j,1} = zeros(numel(idxCluster{j}), 1);
        for k=1:length(idxCluster{j})
            scoreClusters = SegDirScores{j,1}(idxCluster{j,1}{k,1});
            [~, sortIdx] = sort(scoreClusters, 'descend');
            candidate = idxCluster{j,1}{k,1}(sortIdx(1));
            if all(bsxfun(@gt, borderPositions{j}(candidate,:), bboxSmall(:, 1)'),2)...
                    & all(bsxfun(@lt, borderPositions{j}(candidate,:), bboxSmall(:, 2)'),2)
                candidateUse{j,1}(end+1,1) = candidate;
                candidateOutside{j,1}(end+1,1) = 1;
            end
        end
    end
    
    endingsInside = cat(1,candidateOutside{:});
    endingsInside = endingsInside > 0;
    % Logical matrix which indicates endings beyond cutoff border
    endingsOutside = ~endingsInside;
    display([num2str(numel(endingsOutside)) ' endings in total']);
    display([num2str(sum(endingsOutside)) ' endings outside']);
    display([num2str(sum(endingsInside)) ' endings to query']);
    
    % Write out scores and pca, apply all masks
    outside = cellfun('isempty',candidateUse);
    dendrites = dendrites(idxCanidateFound);
    dendrites = dendrites(~redundant);
    dendrites = dendrites(~outside);

    pcaFound = directionality.pca(idxCanidateFound);
    pcaFound = pcaFound(~redundant);
    pcaFound = pcaFound(~outside);
    scoresFound = directionality.scores(idxCanidateFound);
    scoresFound = scoresFound(~redundant);
    scoresFound = scoresFound(~outside);
    borderIdxFound = directionality.borderIdx(idxCanidateFound);
    borderIdxFound = borderIdxFound(~redundant);
    borderIdxFound = borderIdxFound(~outside);
    candidateUse = candidateUse(~outside);
    borderPositions = borderPositions(~outside);

    % Extract some values for borders we want to query
    thisborderPositions = {};
    directions = {};
    for j=1:sum(~outside)
        thisPca{j,1} = [];
        thisScores{j,1} = [];
        thisBorderIdx{j,1} = [];
        thisborderPositions{j,1} = {};
        directions{j,1} = {};
        for k=1:length(candidateUse{j})
            thisPca{j,1}(end+1,:) = pcaFound{j}(:,1,candidateUse{j}(k))';
            thisScores{j,1}(end+1,1) = scoresFound{j}(candidateUse{j}(k));
            thisBorderIdx{j,1}(end+1,1) = borderIdxFound{j}(candidateUse{j}(k));
            thisborderPositions{j,1}{end+1,1} = double(borderPositions{j,1}(candidateUse{j}(k),:));
            directions{j,1}{end+1,1} = bsxfun(@times, bsxfun(@times, sign(thisScores{j,1}(end,1)), thisPca{j,1}(end,:)), 1 ./ param.raw.voxelSize);
        end
    end

    outputFolder = fullfile(dataDir, strcat('dendriteQueries_',num2str(queryRun)));
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder)
    end

    batchBoundaries = round(linspace(1, numel(dendrites), 101));
    for i=1:length(batchBoundaries)-1
        inputList{i,1} = {i};
    end

    sharedInputs = {batchBoundaries,dendrites,thisborderPositions,directions,param,outputFolder};
    sharedInputsLocation = 2:7;
    cluster = Cluster.getCluster('-p -100','-tc 20','-l h_vmem=64G','-l s_rt=5:28:00','-l h_rt=5:29:00');
    job = Cluster.startJob(@connectEM.processQueryTasks, inputList, 'name', 'queryGeneration', 'cluster', cluster, ...
        'sharedInputs', sharedInputs, 'sharedInputsLocation', sharedInputsLocation);

end

