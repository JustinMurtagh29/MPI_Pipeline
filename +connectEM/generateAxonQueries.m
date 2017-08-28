function generateAxonQueries(param)
    % Written by 
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    
    dataDir = fullfile(param.saveFolder, 'aggloState');

    % Load data from ending generation
    endingData = fullfile(dataDir, 'axonEndings.mat');
    endingData = load(endingData);
    idxAll = endingData.out.idxAll;
    idxCanidateFound = endingData.out.axonMask;
    borderIds = endingData.out.borderIds;
    Clusters = endingData.out.borderClusters ;
    
    % Load directionality information
    directionality = fullfile(dataDir, 'axonEndingInputData.mat');
    directionality = load(directionality, 'directionality');
    directionality = directionality.directionality;
    
    % Load border CoMs
    borderCoM = fullfile(dataDir, 'globalBorder.mat');
    borderCoM = load(borderCoM, 'borderCoM');
    borderCoM = borderCoM.borderCoM;
    borderPositions = cellfun(@(x) borderMeta.borderCoM(x,:), directionality.borderIdx(idxCanidateFound),'uni',0);

    % Determine endings which are not redundant(already attached by flight path)
    out = load(fullfile(dataDir, 'axonEndingOverlaps.mat'));
    startEndings = unique(cell2mat(out.startEndingOverlaps));
    endEndings = unique(cell2mat(out.endEndingOverlaps));
    attachedEndings = union(startEndings,endEndings);

    % Load larger 5 micron agglomerates
    m = load(fullfile(dataDir, 'axons_04.mat'));
    axons = m.axons(m.indBigAxons);
    axons = arrayfun(@Agglo.fromSuperAgglo, axons, 'UniformOutput', false);
    clear m
    
    % Write new segmentation based on axon queries
    mapping = connectEM.createLookup(segmentMeta, axonsNew);
    Seg.Global.applyMappingToSegmentation(p, mapping, [outputFolder '1/']);
    clear mapping;
    
    % Find indices of ending candidates in directionality lists (scores,
    % pca...) and exclude redundant endings
    idxUse = idxAll(idxCanidateFound);
    counter = 1;
    idxCluster = [];
    for j=1:length(Clusters)
        idxCluster{j,1} = [];
        for k=1:max(Clusters{j})
            if ~ismember(counter,attachedEndings)
                idxCluster{j,1}{k,1} = [];
            else
                idxCluster{j,1}{k,1} = idxUse{j,1}(find(Clusters{j}==k));
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
    SegDirScores = SegDirScores(redundant);
    borderPositions = borderPositions(redundant);
    
    % Dataset border conditions
    borderNm = repmat(options.border, 1, 3);
    borderVoxel = round(bsxfun(@times, 1./p.raw.voxelSize, borderNm));
    bboxSmall = p.bbox + borderVoxel';

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
            end
        end
    end
    
    % Write out scores and pca, apply all masks
    outside = cellfun('isempty',candidateUse);
    axons = axons(idxCanidateFound);
    axons = axons(redundant);
    axons = axons(~outside);
    
    pcaFound = directionality.pca(idxCanidateFound);
    pcaFound = pcaFound(redundant);
    pcaFound = pcaFound(~outside);
    scoresFound = directionality.scores(idxCanidateFound);
    scoresFound = scoresFound(redundant);
    scoresFound = scoresFound(~outside);
    borderIdxFound = directionality.borderIdx(idxCanidateFound);
    borderIdxFound = borderIdxFound(redundant);
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
            directions{j,1}{end+1,1} = bsxfun(@times, bsxfun(@times, sign(thisScores{j,1}(end,1)), thisPca{j,1}(end,:)), 1 ./ p.raw.voxelSize);
        end
    end
    
    outputFolder = fullfile(param.saveFolder, 'aggloState/');
        
    batchBoundaries = round(linspace(1, numel(axons), 101));
    for i=1:length(batchBoundaries)-1
        inputList{i,1} = {i,batchBoundaries,axons,thisborderPositions,directions,p,outputFolder};
    end
    
    cluster = Cluster.getCluster('-p -100','-tc 20','-l h_vmem=64G','-l s_rt=5:28:00','-l h_rt=5:29:00');
    job = Cluster.startJob(@connectEM.processQueryTasks, inputList, 'name','queryGeneration','cluster', cluster)

    %{
    for i=1:length(batchBoundaries)-1
        tic;
        theseIdx = batchBoundaries(i):batchBoundaries(i+1)-1;
        theseAxons = axons(theseIdx);
        theseBorderPositions = thisborderPositions(theseIdx);
        theseDirections = directions(theseIdx);
        
        q.pos = {};
        q.dir = {};
        q.angles = {};
        idxDelete = [];
        for j=1:length(theseAxons)
            tic;
            thesePositions = cellfun(@(x,y,z)connectEM.correctQueryLocationToEndOfSegment(p, x, y, z, 200), ...
                repmat(theseAxons(j),size(theseBorderPositions{j},1),1), theseBorderPositions{j}, theseDirections{j}, 'uni', 0);
            
            q.pos{end+1} = thesePositions;
            q.dir{end+1} = theseDirections{j}(~cellfun('isempty',thesePositions));
            [phi, theta, psi] = cellfun(@(x)connectEM.calculateEulerAngles(x, p.raw.voxelSize), q.dir{end});
            q.angles{end+1} = mat2cell(cat(2, phi, theta, psi), ones(numel(phi),1), 3);
            toc;
        end
        
        for j=1:length(idxDelete)
            theseAxons{idxDelete(j)} = [];
        end
        theseAxons = theseAxons(~logical(cellfun('isempty',theseAxons)));

        % Calculate euler angles & put into old format
        save([outputFolder 'batch' num2str(i, '%.4i') '.mat'], 'q', 'theseAxons');
        display(['Batch ' num2str(i, '%.4i') ' done']);
        clear these* q phi theta psi;
        toc;
    end
    %}