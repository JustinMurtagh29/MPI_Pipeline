function generateQueriesOfBorderWholeCells(param,suffix,graphInput,runID)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>

    if ~exist('suffix','var')
        suffix = '';
    end

    dataDir = fullfile(param.saveFolder, 'aggloState');

    if nargin < 3
        [graph, ~, borderMeta, ~] = ...
            connectEM.loadAllSegmentationData(param);
    else
        graph = graphInput.graph;
        borderMeta = graphInput.borderMeta;
    end
   
    % Load data from ending generation
    endingData = fullfile(dataDir, sprintf('borderWholeCellsEndingsAllData_%s.mat',suffix));
    endingData = load(endingData);
    % Extract some variables
    idxAll = endingData.idxAll;
    idxCanidateFound = endingData.candidateMask;
    borderIds = endingData.borderIds;
    borderClusters = endingData.borderClusters ;

    % Load directionality information
    directionality = fullfile(dataDir, sprintf('borderWholeCellsEndingInputData_%s.mat',suffix));
    directionality = load(directionality, 'directionality');
    directionality = directionality.directionality;
    
    % Load border CoMs
    borderCoM = borderMeta.borderCoM;
    borderPositions = cellfun(@(x) borderCoM(x,:), directionality.borderIdx(idxCanidateFound),'uni',0);

    % Load larger 5 micron agglomerates
    m = load(fullfile(dataDir, 'dendrites_08.mat'));
    borderWholeCells = m.dendrites(m.BorderWholeCellId);
    superDendrites = borderWholeCells;
    wholeCells = arrayfun(@Agglo.fromSuperAgglo, borderWholeCells, 'UniformOutput', false);
    
    % Find indices of ending candidates in directionality lists (scores,
    % pca...)
    idxUse = idxAll(idxCanidateFound);
    counter = 1;
    idxCluster = [];
    for j=1:length(borderClusters)
        idxCluster{j,1} = {};
        for k=1:max(borderClusters{j})
            idxCluster{j,1}{k,1} = idxUse{j,1}(find(borderClusters{j}==k));
            counter = counter + 1;
        end
        idxCluster{j,1} = idxCluster{j,1}(~cellfun('isempty',idxCluster{j,1}));
    end

    % Write out absolut values of seg direction scores
    SegDirScores = cellfun(@(x)abs(x),directionality.scores,'uni',0);
    
    % Dataset border conditions
    options.border = [3000; -3000];
    borderNm = repmat(options.border, 1, 3);
    borderVoxel = round(bsxfun(@times, 1./param.raw.voxelSize, borderNm));
    bboxSmall = param.bbox + borderVoxel';

    % Keep only candidate with highest score for each cluster and exclude
    % candidate outside border cutoff
    candidateUse = {};
    candidateInside = {};
    for j=1:length(idxCluster)
        candidateUse{j,1} = [];
        candidateInside{j,1} = zeros(numel(idxCluster{j}), 1);
        for k=1:length(idxCluster{j})
            scoreClusters = SegDirScores{j,1}(idxCluster{j,1}{k,1});
            [~, sortIdx] = sort(scoreClusters, 'descend');
            candidate = idxCluster{j,1}{k,1}(sortIdx(1));
            if all(bsxfun(@gt, borderPositions{j}(candidate,:), bboxSmall(:, 1)'),2)...
                    & all(bsxfun(@lt, borderPositions{j}(candidate,:), bboxSmall(:, 2)'),2)
                candidateUse{j,1}(end+1,1) = candidate;
                candidateInside{j,1}(k,1) = 1;
            end
        end
    end
    
    endingsInside = cat(1,candidateInside{:});
    endingsInside = endingsInside > 0;
    % Logical matrix which indicates endings beyond cutoff border
    endingsOutside = ~endingsInside;
    display([num2str(numel(endingsOutside)) ' endings in total']);
    display([num2str(sum(endingsOutside)) ' endings outside']);
    display([num2str(sum(endingsInside)) ' endings to query']);
    
    % Write out scores and pca, apply all masks
    outside = cellfun('isempty',candidateUse);
    wholeCells = wholeCells(~outside);

    pcaFound = directionality.pca;
    pcaFound = pcaFound(~outside);
    scoresFound = directionality.scores;
    scoresFound = scoresFound(~outside);
    borderIdxFound = directionality.borderIdx;
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

    borderEdges = graph.edges(~isnan(graph.borderIdx), :);
    
    outputFolder = fullfile(dataDir, strcat('borderWholeCellDendriteQueries/wholeCellDendriteQueries_',num2str(runID),'/'));
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder)
    end
    
    superDendrites = superDendrites(~outside);
    IDs = find(~outside);
    
    parameters.experiment.name='2012-09-28_ex145_07x2_ROI2017_dendrites_20171018';
    parameters.scale.x = '11.24';
    parameters.scale.y = '11.24';
    parameters.scale.z = '28';
    parameters.offset.x = '0';
    parameters.offset.y = '0';
    parameters.offset.z = '0';
    
    for i=1:length(IDs)
        theseBorderEdges =  borderEdges(thisBorderIdx{i},:);
        aggloEndingSegIds = theseBorderEdges(ismember(theseBorderEdges,wholeCells{i}));
        
        % Filtering of endings with a neighboring flight path or mistakenly
        % detected because of missing segment.
        % Edges of superagglo in SegIDs
        segEdges = zeros(size(superDendrites(i).edges));
        segEdges(:,1) = superDendrites(i).nodes(superDendrites(i).edges(:,1),4);
        segEdges(:,2) = superDendrites(i).nodes(superDendrites(i).edges(:,2),4);
        % Neighbouring edges of ending segments
        [segLocRow,~] = arrayfun(@(x)find(ismember(segEdges,x)),aggloEndingSegIds,'uni',0);
        neighbours = cellfun(@(x)unique(segEdges(x,:)),segLocRow,'uni',0);
        for j=1:length(aggloEndingSegIds)
            neighbours{j}(neighbours{j} == aggloEndingSegIds(j)) = [];
        end
        % All surrounded segments of ending segment.
        allNeighbours = arrayfun(@(x)graph.neighbours(x),aggloEndingSegIds);
        surroundingCheck = cellfun(@(x,y)ismember(x,y),neighbours,allNeighbours,'uni',0);
        
        flightPathCheck = cellfun(@(x)any(isnan(x)),neighbours);
        skippedSegmentCheck = cellfun(@(x)any(x==0),surroundingCheck);
        aggloEndingSegIds(flightPathCheck|skippedSegmentCheck) = [];
        
        idxComments = ismember(superDendrites(i).nodes(:,4),aggloEndingSegIds);
        superDendrites(i).comments = repmat({''},size(superDendrites(i).nodes,1),1);
        superDendrites(i).comments(idxComments) = repmat({'ending'},sum(idxComments),1);
    
        if ~any(idxComments)
            continue
        end
        
        connectEM.generateSkeletonFromAggloNew(superDendrites(i), {sprintf('wholeCellAgglo_%02d',i)} , outputFolder, [],parameters,sprintf('WholeCell%s_%02d.nml',suffix,IDs(i)));
    end

end

