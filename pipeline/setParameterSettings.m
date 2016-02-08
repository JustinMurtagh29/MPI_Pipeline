function [p, pT] = setParameterSettings(p)
    % Pass structure with basic settings to add all dependent and constant p for pipeline

    % Size of local segmentation and local graph construction
    p.tileSize =  [512; 512; 256];
    % Check whether bounding box meets requirements and fix if not
    p.bbox = fixBoundingBox(p);
    % Add a temp folder for storing intermediate results
    p.tempFolder = [p.saveFolder 'temp/'];
    % Overlap between tiles in segmentation for gloablization
    p.tileBorder = [-256 256; -256 256; -128 128];
    p.tiles = (p.bbox(:,2) - p.bbox(:,1) + 1) ./ p.tileSize;
    % Which function to use to normalize data to zero mean and one std
    [meanVal, stdVal] = determineMeanAndStdOfData(p);
    p.norm.func = @(x)normalizeStack(x,meanVal,stdVal);
    % Which classifier to use
    p.cnn.dateStrings = '20130516T204040';
    p.cnn.iter = 8; 
    p.cnn.gpu = 3;
    p.cnn.first = ['/gaba/u/mberning/results/parameterSearch/' p.cnn.dateStrings ...
        '/iter' num2str(p.cnn.iter, '%.2i') '/gpu' num2str(p.cnn.gpu, '%.2i') '/saveNet0000000001.mat'];
    p.cnn.GPU = false;
    % Function to use for classification
    p.class.func = @bigFwdPass;
    % Location to store CNN classification
    p.class.root = [p.saveFolder 'class/'];
    p.class.prefix = p.raw.prefix;
    % Function to use for segmentation
    p.seg.func = @(x)watershedSeg_v1_cortex(x,{p.seg.threshold 10});
    p.seg.root = [p.saveFolder 'globalSeg/'];
    p.seg.prefix = p.raw.prefix;
    % Specify arguments for filterbank applied to raw and aff data each
    p.filter = {{'sortedeigenvalueshessian' [3 5] []}...
        {'gaussiansmoothedgradmagnitude' [3 5] []}...
        {'intensitygaussiansmoothed' [3 5] []}...
        {'sortedeigenvaluesstructure' [3 5] [5 7]}...
        {'laplaceofgaussian' [3 5] []}...
        {'differenceofgaussians' [3 5] []}};
    % Feature p
    p.feature.root = [p.saveFolder 'features/'];
    % Function to use for feature calculation
    p.feature.func = @calcFeatures;
    % Choose to filter 'raw' and 'class' data
    p.feature.input = {'raw', 'aff'};
    % Correspondence p
    p.correspondence.overlap = 1; % overlap of local segmentation to compare on each side around a face
    p.correspondence.saveFolder = [p.tempFolder 'correspondences/'];

    % GLOBAL SETTINGS FOR fromGraphToDB.m
    % State variables from the GP
    p.gp.stateFolder = [p.saveFolder 'state/'];
    p.gp.normValues = [p.gp.stateFolder 'normValues.mat'];
    p.gp.hyperParameter = [p.gp.stateFolder 'hyperParameter.mat'];
    p.gp.initalGroundTruth = [p.gp.stateFolder 'initalGroundTruth.mat'];
    % Define cutoff(s) for writing to knowledge DB 
    p.gp.upperCut = .95;
    p.gp.lowerCut = .15;

    % LOCAL SETTINGS for each tile
    for i=1:p.tiles(1)
        for j=1:p.tiles(2)
            for k=1:p.tiles(3)
                % Save path for data relating to this tile
                p.local(i,j,k).saveFolder = [p.saveFolder 'local/' ...
                    'x' num2str(i, '%.4i') 'y' num2str(j, '%.4i') 'z' num2str(k, '%.4i') '/'];
                % Bounding box without border for this tile
                p.local(i,j,k).bboxSmall = [p.bbox(:,1) + [i-1; j-1; k-1] .* p.tileSize p.bbox(:,1) ...
                    + [i; j; k] .* p.tileSize - [1; 1; 1]];
                % Restrict bounding box big in order to avoid overlaps at border cubes of segmentation
                bboxBig =  p.local(i,j,k).bboxSmall + p.tileBorder;
                bboxBig(:,1) = arrayfun(@(x,y)max(x,y),bboxBig(:,1),p.bbox(:,1));
                bboxBig(:,2) = arrayfun(@(x,y)min(x,y),bboxBig(:,2),p.bbox(:,2));
                p.local(i,j,k).bboxBig = bboxBig;
                % Where to save
                p.local(i,j,k).segFile = [p.local(i,j,k).saveFolder 'seg.mat'];
                p.local(i,j,k).tempSegFile = strrep(p.local(i,j,k).segFile, '/local/', '/temp/');
                p.local(i,j,k).edgeFile = [p.local(i,j,k).saveFolder 'edges.mat'];
                p.local(i,j,k).borderFile =  [p.local(i,j,k).saveFolder 'borders.mat'];
                p.local(i,j,k).weightFile = [p.local(i,j,k).saveFolder 'weights.mat'];
                p.local(i,j,k).probFile = [p.local(i,j,k).saveFolder 'prob.mat'];
                % Same files for glia prediction
                p.local(i,j,k).segmentFile = [p.local(i,j,k).saveFolder 'segments.mat'];    
                p.local(i,j,k).segmentWeightFile = [p.local(i,j,k).saveFolder 'segmentWeights.mat'];
                p.local(i,j,k).gliaProbFile = [p.local(i,j,k).saveFolder 'gliaProb.mat'];
            end
        end
    end


    % GLOBAL SETTINGS FOR training data generation
    pT = p;
    % Remove all fields that do not make sense in training data setting
    pT = rmfield(pT, {'local' 'bbox' 'tileSize' 'tiles' 'correspondence'});
    pT.cnn.GPU = false; % minicubeFwdPass on CPU to allow arbitrary region size

    % Densly skeletonized regions in dataset
    % Region from Heiko
    pT.local(1).bboxSmall = [4097 4736; 4481 5248; 2250 2450];
    pT.local(1).trainFileRaw = '/gaba/u/mberning/data/cortex/denseSkel/region1.nml';
    pT.local(1).trainFileGlia = '/gaba/u/mberning/data/cortex/denseSkel/region1glia.nml';
    pT.local(1).trainFileLocal = '/gaba/u/mberning/data/cortex/denseSkel/region1local.nml'; 
    pT.local(1).trainFileLocalWithoutGlia = '/gaba/u/mberning/data/cortex/denseSkel/region1localWithoutGlia.nml';
    % Region from Alex
    pT.local(2).bboxSmall = [1417 1717; 4739 5039; 890 1190];
    pT.local(2).trainFileRaw = '/gaba/u/mberning/data/cortex/denseSkel/region2.nml';
    pT.local(2).trainFileGlia = '/gaba/u/mberning/data/cortex/denseSkel/region2glia.nml';
    pT.local(2).trainFileLocal = '/gaba/u/mberning/data/cortex/denseSkel/region2local.nml'; 
    pT.local(2).trainFileLocalWithoutGlia = '/gaba/u/mberning/data/cortex/denseSkel/region2localWithoutGlia.nml';
    % Region from Max & Anna
    pT.local(3).bboxSmall = [6800 7100; 2140 2440; 1236 1536];
    pT.local(3).trainFileRaw = '/gaba/u/mberning/data/cortex/denseSkel/region3.nml';
    pT.local(3).trainFileGlia = '/gaba/u/mberning/data/cortex/denseSkel/region3glia.nml';
    pT.local(3).trainFileLocal = '/gaba/u/mberning/data/cortex/denseSkel/region3local.nml';
    pT.local(3).trainFileLocalWithoutGlia = '/gaba/u/mberning/data/cortex/denseSkel/region3localWithoutGlia.nml';

    % LOCAL SETTINGS FOR training tiles
    for i=1:3
        % Save path for data relating to this tile
        pT.local(i).saveFolder = [pT.saveFolder 'train' num2str(i, '%.4i') '/'];
        % Bounding box without and with border for this tile
        pT.local(i).bboxBig = pT.local(i).bboxSmall + pT.tileBorder;
        % Where to save
        pT.local(i).class.root = [pT.local(i).saveFolder 'class/'];
        pT.local(i).class.prefix = pT.class.prefix;
        pT.local(i).seg.pSearchFolder = [pT.local(i).saveFolder 'pSearch/'];
        pT.local(i).segFile = [pT.local(i).saveFolder 'seg.mat'];
        pT.local(i).edgeFile = [pT.local(i).saveFolder 'edges.mat'];
        pT.local(i).borderFile =  [pT.local(i).saveFolder 'borders.mat'];
        pT.local(i).weightFile = [pT.local(i).saveFolder 'weights.mat'];
        pT.local(i).gtFile = [pT.local(i).saveFolder 'region' num2str(i) '.mat'];
        % Benjamin's glia prediction
        pT.local(i).segmentFile = [pT.local(i).saveFolder 'segments.mat'];    
        pT.local(i).segmentWeightFile = [pT.local(i).saveFolder 'segmentWeights.mat'];
        pT.local(i).gliaProbFile = [pT.local(i).saveFolder 'gliaProb.mat'];
    end

    % Create state folder for this p settings GP
    if ~exist(p.gp.stateFolder, 'dir')
        mkdir(p.gp.stateFolder);
    end

    % Save everything
    save([p.saveFolder 'allParameter.mat'], 'p', 'pT');

end

function bbox = fixBoundingBox(p)
    % Make bounding box meet requirements

    % Transform from webKNOSSOS style bounding box to bounding box format
    % used in pipeline
    bbox = reshape(p.bbox_wK,3,2);

    % First check whether bounding box is aligned with KNOSSOS cubes
    lowerLimitMod = mod(bbox(:,1)-1,128);
    upperLimitMod = mod(bbox(:,2),128);
    if any(lowerLimitMod)
        warning('Lower edge of bounding box not aligned to KNOSSOS cubes, fixing.');
        bbox(:,1) = bbox(:,1) + mod(128 - lowerLimitMod, 128);
    end
    if any(upperLimitMod)
        warning('Upper edge of bounding box not aligned to KNOSSOS cubes, fixing.');
        bbox(:,2) = bbox(:,2) - upperLimitMod;
    end

    % Second check whether it is 'tileable' using tileSize
    sizeOfRoi = bbox(:,2) - bbox(:,1) + 1;
    tileSize = p.tileSize;
    sizeToRemove = mod(sizeOfRoi, tileSize);
    if any(sizeToRemove)
        warning(['Bounding box not divisable by tileSize (' num2str(tileSize') '), fixing']);
        for i=1:3
            if sizeToRemove(i)
                if ~mod(sizeToRemove(i),2*128)
                    % If sizeToRemove is divisiable by 2 * knossos cube size, remove at both ends
                    bbox(i,1) = bbox(i,1) + sizeToRemove(i)/2;
                    bbox(i,2) = bbox(i,2) - sizeToRemove(i)/2;
                else
                    % Otherwise arbitrarily drop more at upper edge
                    bbox(i,1) = bbox(i,1) + (sizeToRemove(i)-128)/2;
                    bbox(i,2) = bbox(i,2) - ((sizeToRemove(i)-128)/2+128);
                end
            end
        end
    end

end

function [meanVal, stdVal] = determineMeanAndStdOfData(p)

    display('Sampling mean and standard deviation values for CNN normalization');
    % How many 100^3 samples to use for determination of normalization
    % constants, 100 seems rather too much but ok as only takes 5 min
    nrCubesToSample = 100;
    sizeOfRoi = p.bbox(:,2) - p.bbox(:,1) + 1;
    meanVal = zeros(nrCubesToSample,1);
    stdVal = zeros(nrCubesToSample,1);
    for i=1:nrCubesToSample
        lowerLeft = [randi(sizeOfRoi(1)-99); randi(sizeOfRoi(2)-99); randi(sizeOfRoi(3)-99)];
        bbox = cat(2,lowerLeft, lowerLeft + 99);
        raw = loadRawData(p.raw.root, p.raw.prefix, bbox, false);
        meanVal(i) = mean(raw(:));
        stdVal(i) = std(raw(:));
    end

    meanVal = median(meanVal);
    stdVal = median(stdVal);

end
