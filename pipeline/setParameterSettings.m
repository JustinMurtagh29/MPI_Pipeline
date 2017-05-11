function [p, pT] = setParameterSettings(p)
    % Pass structure with basic settings to add all dependent and constant p for pipeline

    % Sanitize paths by adding trailing slashes
    if p.saveFolder(end) ~= filesep
        p.saveFolder = [p.saveFolder, filesep];
        warning('Added trailing slash to p.saveFolder');
    end

    if p.raw.root(end) ~= filesep
        p.raw.root = [p.raw.root, filesep];
        warning('Added trailing slash to p.raw.root');
    end

    % Size of local segmentation and local graph construction
    p.tileSize =  [512; 512; 256];
    
    % Check whether bounding box meets requirements and fix if not
    p.bbox = fixBoundingBox(p);
    
    % Add a temp folder for storing intermediate results
    p.tempFolder = [p.saveFolder 'temp/'];
    
    % Overlap between tiles in segmentation for gloablization
    p.tileBorder = [-256 256; -256 256; -128 128];
    p.tiles = ceil((p.bbox(:,2) - p.bbox(:,1) + 1) ./ p.tileSize);
    
    % Which function to use to normalize data to zero mean and one std
    [meanVal, stdVal] = ...
        Knossos.estGlobalMeanAndStd(p);
    p.norm.func = @(x)normalizeStack(x,meanVal,stdVal);
    
    % Which classifier to use
    p.cnn.dateStrings = '20130516T204040';
    p.cnn.iter = 8; 
    p.cnn.gpu = 3;
    p.cnn.first = ['/gaba/u/mberning/results/parameterSearch/' p.cnn.dateStrings ...
        '/iter' num2str(p.cnn.iter, '%.2i') '/gpu' num2str(p.cnn.gpu, '%.2i') '/saveNet0000000001.mat'];
    p.cnn.GPU = false;
    
    % Location to store CNN classification
    p.class = Util.modifyStruct( ...
        p.raw, ...
        'func', @bigFwdPass, ...
        'root', strcat(p.tempFolder, 'class', filesep));
    
    % Function to use for segmentation
    p.seg = Util.modifyStruct( ...
        p.raw, ...
        'func', @(x) watershedSeg_v1_cortex(x, {p.seg.threshold, 10}), ...
        'root', strcat(p.saveFolder, 'globalSeg', filesep), ...
        'threshold', p.seg.threshold);
    
    % Specify arguments for filterbank applied to raw and aff data each
    p.filter = {
        {'sortedeigenvalueshessian' [3 5] []}...
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
    p.gp.initialGroundTruth = [p.gp.stateFolder 'initialGroundTruth.mat'];
    
    % Define cutoff(s) for writing to knowledge DB 
    %p.gp.upperCut = .95;
    %p.gp.lowerCut = .15;

    % LOCAL SETTINGS for each tile
    for i=1:p.tiles(1)
        for j=1:p.tiles(2)
            for k=1:p.tiles(3)
                % Save path for data relating to this tile
                p.local(i,j,k).saveFolder = [p.saveFolder 'local/' ...
                    'x' num2str(i, '%.4i') 'y' num2str(j, '%.4i') 'z' num2str(k, '%.4i') '/'];
                
                % Bounding box without border for this tile
                p.local(i,j,k).bboxSmall = [ ...
                    p.bbox(:,1) + [i-1; j-1; k-1] .* p.tileSize, ...
                    p.bbox(:,1) + [i;   j;   k  ] .* p.tileSize - [1; 1; 1]];
                p.local(i,j,k).bboxSmall(:, 2) = ...
                    min(p.local(i,j,k).bboxSmall(:, 2), p.bbox(:, 2));
                                        
                % Restrict bounding box big in order to avoid overlaps at border cubes of segmentation
                bboxBig =  p.local(i,j,k).bboxSmall + p.tileBorder;
                bboxBig(:,1) = max(bboxBig(:,1),p.bbox(:,1));
                bboxBig(:,2) = min(bboxBig(:,2),p.bbox(:,2));
                p.local(i,j,k).bboxBig = bboxBig;
                
                % Where to save
                p.local(i,j,k).segFile = [p.local(i,j,k).saveFolder 'seg.mat'];
                p.local(i,j,k).tempSegFile = strrep(p.local(i,j,k).segFile, '/local/', '/temp/');
                p.local(i,j,k).edgeFile = [p.local(i,j,k).saveFolder 'edges.mat'];
                p.local(i,j,k).borderFile =  [p.local(i,j,k).saveFolder 'borders.mat'];
                p.local(i,j,k).weightFile = [p.local(i,j,k).saveFolder 'weights.mat'];
                p.local(i,j,k).probFile = [p.local(i,j,k).saveFolder 'prob.mat'];
                p.local(i,j,k).synapseFile = [p.local(i,j,k).saveFolder 'synapses.mat'];

                % Same files for glia prediction
                p.local(i,j,k).segmentFile = [p.local(i,j,k).saveFolder 'segments.mat'];    
                p.local(i,j,k).segmentWeightFile = [p.local(i,j,k).saveFolder 'segmentWeights.mat'];
                p.local(i,j,k).gliaProbFile = [p.local(i,j,k).saveFolder 'gliaProb.mat'];
            end
        end
    end

    % Create state folder for this p settings GP
    if ~exist(p.gp.stateFolder, 'dir')
        mkdir(p.gp.stateFolder);
    end

    % Save everything
    Util.save([p.saveFolder 'allParameter.mat'], p);

end

function bbox = fixBoundingBox(p)
    % Make bounding box meet requirements

    % Transform from webKNOSSOS style bounding box to
    % bounding box format used in pipeline
    bbox = reshape(p.bbox_wK, 3, 2);
    % webKnossos now uses width height depth 
    bbox(:,1) = bbox(:,1) + 1;
    bbox(:,2) = bbox(:,1) + bbox(:,2) - 1;

    % First check whether bounding box is aligned with KNOSSOS cubes
    lowerLimitMod = mod(bbox(:,1) - 1, 128);
    upperLimitMod = mod(diff(bbox, 1, 2) + 1, p.tileSize);


    
    if any(lowerLimitMod)
        warning('Lower edge of bounding box not aligned to KNOSSOS cubes, fixing, please check p.bbox');
        idx = lowerLimitMod ~= 0;
        bbox(idx,1) = bbox(idx,1) + 128 - lowerLimitMod;
    end
    
    if any(upperLimitMod < 64 & upperLimitMod ~= 0)
        error('Upper edge of bounding box produces small last cube.');
    end
end

