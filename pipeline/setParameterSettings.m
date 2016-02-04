function [parameter, parameterTrain] = setParameterSettings(parameter)
    % Pass structure with basic settings to add all dependent and constant parameter for pipeline

    % Size of local segmentation and local graph construction
    parameter.tileSize =  [512; 512; 256];
    % Check whether bounding box meets requirements and fix if not
    parameter.bbox = fixBoundingBox(parameter);
    % Overlap between tiles in segmentation for gloablization
    parameter.tileBorder = [-256 256; -256 256; -128 128];
    parameter.tiles = (parameter.bbox(:,2) - parameter.bbox(:,1) + 1) ./ parameter.tileSize;
    % Which function to use to normalize data to zero mean and one std
    [meanVal, stdVal] = determineMeanAndStdOfData(parameter);
    parameter.norm.func = @(x)normalizeStack(x,meanVal,stdVal);
    % Which classifier to use
    parameter.cnn.dateStrings = '20130516T204040';
    parameter.cnn.iter = 8; 
    parameter.cnn.gpu = 3;
    parameter.cnn.first = ['/gaba/u/mberning/results/parameterSearch/' parameter.cnn.dateStrings ...
        '/iter' num2str(parameter.cnn.iter, '%.2i') '/gpu' num2str(parameter.cnn.gpu, '%.2i') '/saveNet0000000001.mat'];
    parameter.cnn.GPU = true;
    % Function to use for classification
    parameter.class.func = @bigFwdPass;
    % Location to store CNN classification
    parameter.class.root = [parameter.saveFolder 'class/'];
    parameter.class.prefix = parameter.raw.prefix;
    % Function to use for segmentation
    parameter.seg.func = @(x)watershedSeg_v1_cortex(x,{parameter.seg.threshold 10});
    parameter.seg.root = [parameter.saveFolder 'globalSeg/'];
    parameter.seg.prefix = parameter.raw.prefix;
    % Specify arguments for filterbank applied to raw and aff data each
    parameter.filter = {{'sortedeigenvalueshessian' [3 5] []}...
        {'gaussiansmoothedgradmagnitude' [3 5] []}...
        {'intensitygaussiansmoothed' [3 5] []}...
        {'sortedeigenvaluesstructure' [3 5] [5 7]}...
        {'laplaceofgaussian' [3 5] []}...
        {'differenceofgaussians' [3 5] []}};
    % Feature parameter
    parameter.feature.root = [parameter.saveFolder 'features/'];
    % Function to use for feature calculation
    parameter.feature.func = @calcFeatures;
    % Choose to filter 'raw' and 'class' data
    parameter.feature.input = {'raw', 'aff'};
    % Correspondence parameter
    parameter.correspondence.overlap = 1; % overlap of local segmentation to compare on each side around a face
    parameter.correspondence.saveFolder = [parameter.saveFolder 'correspondences/'];

    % GLOBAL SETTINGS FOR fromGraphToDB.m
    % State variables from the GP
    parameter.gp.stateFolder = [parameter.saveFolder 'state/'];
    parameter.gp.normValues = [parameter.gp.stateFolder 'normValues.mat'];
    parameter.gp.hyperParameter = [parameter.gp.stateFolder 'hyperParameter.mat'];
    parameter.gp.initalGroundTruth = [parameter.gp.stateFolder 'initalGroundTruth.mat'];
    % Define cutoff(s) for writing to knowledge DB 
    parameter.gp.upperCut = .95;
    parameter.gp.lowerCut = .15;

    % LOCAL SETTINGS for each tile
    for i=1:parameter.tiles(1)
        for j=1:parameter.tiles(2)
            for k=1:parameter.tiles(3)
                % Save path for data relating to this tile
                parameter.local(i,j,k).saveFolder = [parameter.saveFolder 'local/' ...
                    'x' num2str(i, '%.4i') 'y' num2str(j, '%.4i') 'z' num2str(k, '%.4i') '/'];
                % Bounding box without and with border for this tile
                parameter.local(i,j,k).bboxSmall = [parameter.bbox(:,1) + [i-1; j-1; k-1] .* parameter.tileSize parameter.bbox(:,1) ...
                    + [i; j; k] .* parameter.tileSize - [1; 1; 1]];
                parameter.local(i,j,k).bboxBig = parameter.local(i,j,k).bboxSmall + parameter.tileBorder;
                % Where to save
                parameter.local(i,j,k).segFile = [parameter.local(i,j,k).saveFolder 'seg.mat'];
                parameter.local(i,j,k).edgeFile = [parameter.local(i,j,k).saveFolder 'edges.mat'];
                parameter.local(i,j,k).borderFile =  [parameter.local(i,j,k).saveFolder 'borders.mat'];
                parameter.local(i,j,k).weightFile = [parameter.local(i,j,k).saveFolder 'weights.mat'];
                parameter.local(i,j,k).probFile = [parameter.local(i,j,k).saveFolder 'prob.mat'];
                % Same files for glia prediction
                parameter.local(i,j,k).segmentFile = [parameter.local(i,j,k).saveFolder 'segments.mat'];    
                parameter.local(i,j,k).segmentWeightFile = [parameter.local(i,j,k).saveFolder 'segmentWeights.mat'];
                parameter.local(i,j,k).gliaProbFile = [parameter.local(i,j,k).saveFolder 'gliaProb.mat'];
            end
        end
    end


    % GLOBAL SETTINGS FOR training data generation
    parameterTrain = parameter;
    % Remove all fields that do not make sense in training data setting
    parameterTrain = rmfield(parameterTrain, {'local' 'bbox' 'tileSize' 'tiles' 'correspondence'});
    parameterTrain.cnn.GPU = false; % minicubeFwdPass on CPU to allow arbitrary region size

    % Densly skeletonized regions in dataset
    % Region from Heiko
    parameterTrain.local(1).bboxSmall = [4097 4736; 4481 5248; 2250 2450];
    parameterTrain.local(1).trainFileRaw = '/gaba/u/mberning/data/cortex/denseSkel/region1.nml';
    parameterTrain.local(1).trainFileGlia = '/gaba/u/mberning/data/cortex/denseSkel/region1glia.nml';
    parameterTrain.local(1).trainFileLocal = '/gaba/u/mberning/data/cortex/denseSkel/region1local.nml'; 
    parameterTrain.local(1).trainFileLocalWithoutGlia = '/gaba/u/mberning/data/cortex/denseSkel/region1localWithoutGlia.nml';
    % Region from Alex
    parameterTrain.local(2).bboxSmall = [1417 1717; 4739 5039; 890 1190];
    parameterTrain.local(2).trainFileRaw = '/gaba/u/mberning/data/cortex/denseSkel/region2.nml';
    parameterTrain.local(2).trainFileGlia = '/gaba/u/mberning/data/cortex/denseSkel/region2glia.nml';
    parameterTrain.local(2).trainFileLocal = '/gaba/u/mberning/data/cortex/denseSkel/region2local.nml'; 
    parameterTrain.local(2).trainFileLocalWithoutGlia = '/gaba/u/mberning/data/cortex/denseSkel/region2localWithoutGlia.nml';
    % Region from Max & Anna
    parameterTrain.local(3).bboxSmall = [6800 7100; 2140 2440; 1236 1536];
    parameterTrain.local(3).trainFileRaw = '/gaba/u/mberning/data/cortex/denseSkel/region3.nml';
    parameterTrain.local(3).trainFileGlia = '/gaba/u/mberning/data/cortex/denseSkel/region3glia.nml';
    parameterTrain.local(3).trainFileLocal = '/gaba/u/mberning/data/cortex/denseSkel/region3local.nml';
    parameterTrain.local(3).trainFileLocalWithoutGlia = '/gaba/u/mberning/data/cortex/denseSkel/region3localWithoutGlia.nml';

    % LOCAL SETTINGS FOR training tiles
    for i=1:3
        % Save path for data relating to this tile
        parameterTrain.local(i).saveFolder = [parameterTrain.saveFolder 'train' num2str(i, '%.4i') '/'];
        % Bounding box without and with border for this tile
        parameterTrain.local(i).bboxBig = parameterTrain.local(i).bboxSmall + parameterTrain.tileBorder;
        % Where to save
        parameterTrain.local(i).class.root = [parameterTrain.local(i).saveFolder 'class/'];
        parameterTrain.local(i).class.prefix = parameterTrain.class.prefix;
        parameterTrain.local(i).seg.parameterSearchFolder = [parameterTrain.local(i).saveFolder 'parameterSearch/'];
        parameterTrain.local(i).segFile = [parameterTrain.local(i).saveFolder 'seg.mat'];
        parameterTrain.local(i).edgeFile = [parameterTrain.local(i).saveFolder 'edges.mat'];
        parameterTrain.local(i).borderFile =  [parameterTrain.local(i).saveFolder 'borders.mat'];
        parameterTrain.local(i).weightFile = [parameterTrain.local(i).saveFolder 'weights.mat'];
        parameterTrain.local(i).gtFile = [parameterTrain.local(i).saveFolder 'region' num2str(i) '.mat'];
        % Benjamin's glia prediction
        parameterTrain.local(i).segmentFile = [parameterTrain.local(i).saveFolder 'segments.mat'];    
        parameterTrain.local(i).segmentWeightFile = [parameterTrain.local(i).saveFolder 'segmentWeights.mat'];
        parameterTrain.local(i).gliaProbFile = [parameterTrain.local(i).saveFolder 'gliaProb.mat'];
    end

    % Create state folder for this parameter settings GP
    if ~exist(parameter.gp.stateFolder, 'dir')
        mkdir(parameter.gp.stateFolder);
    end

    % Save everything
    pT = parameterTrain;
    p = parameter;
    save([parameter.saveFolder 'allParameter.mat'], 'p', 'pT');

end

function bbox = fixBoundingBox(parameter)
    % Make bounding box meet requirements
    
    % Transform from webKNOSSOS style bounding box to bounding box format
    % used in pipeline
    bbox = reshape(parameter.bbox_wK,3,2);
    
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
    tileSize = parameter.tileSize;
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

function [meanVal, stdVal] = determineMeanAndStdOfData(parameter)

display('Sampling mean and standard deviation values for CNN normalization');
% How many 100^3 samples to use for determination of normalization
% constants, 100 seems rather too much but ok as only takes 5 min
nrCubesToSample = 100;
sizeOfRoi = parameter.bbox(:,2) - parameter.bbox(:,1) + 1;
meanVal = zeros(nrCubesToSample,1);
stdVal = zeros(nrCubesToSample,1);
for i=1:nrCubesToSample
    lowerLeft = [randi(sizeOfRoi(1)-99); randi(sizeOfRoi(2)-99); randi(sizeOfRoi(3)-99)];
    bbox = cat(2,lowerLeft, lowerLeft + 99);
    raw = loadRawData(parameter.raw.root, parameter.raw.prefix, bbox, false);
    meanVal(i) = mean(raw(:));
    stdVal(i) = std(raw(:));
end

meanVal = median(meanVal);
stdVal = median(stdVal);

end
