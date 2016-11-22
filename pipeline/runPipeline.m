function runPipeline(p)

    % Runs CNN based forward pass for region defined in p.bbox,p.raw and saves as p.class	
    % Because CNN is translation invariant, saved as KNOSSOS hierachy again 
    % Pass bounding box as well as tileSize will be added in bigFwdPass otherwise (legacy stuff)
    % This uses CNN subfolder in code repository
    %job = bigFwdPass(p, p.bbox);
    %Cluster.waitForJob(job);
    
    % If you have a lot of myelin in your dataset, you might
    % want to consider applying the following "myelin fix".
    % This is an ugly hack and was originally developed for
    % Jakob's mouse M1 L6 dataset. It will need some tweaking
    % before it works on other datasets. For more info, check
    % out the documentation of +Myelin/runFixInBox
    
    % % define new classification prefix
    % newPrefix = [p.class.prefix, '_mod'];
    %
    % % run myelin detection
    % job = Myelin.runFix(p, newPrefix);
    % Cluster.waitForJob(job);
    %
    % % use modified classification in subsequent steps
    % p.class.prefix = newPrefix;
    % clear newPrefix;
    
    % Runs watershed based segmentation for region defined in p.bbox on p.raw and saves as p.local.segFile
    % Because watershed segmentation has FOV effects (no translation invariance), processed with large
    % overlap and later joined together (see correspondences)
    % Uses segmentation subfolder in code repository
    %job = miniSegmentation(p);
    %Cluster.waitForJob(job);
    
    % Find correspondences between tiles processed (segmented) in paralell
    % Uses correspondences subfolder in code repository
    %job = correspondenceFinder(p);
    %Cluster.waitForJob(job);
    
    % Transfer segmentation from tempSegFile to segFile and drop overlaps 
    % Also in correspondences subfolder
    %job = removeOverlaps(p);
    %Cluster.waitForJob(job);
    
    % Make segmentation IDs unique over dataset (were only unique in each
    % tile before), this will be called global IDs from time to time
    % will work on correspondences and segmentation
    % see globalization subfolder 
    % globalization is run in two parts: globalizeSegmentation and globalizeCorrespondences
    %job = globalizeSegmentation(p);
    %Cluster.waitForJob(job);
    
    %job = globalizeCorrespondences(p);
    %Cluster.waitForJob(job);
    
    % Construct graph on globalized version of segmentation
    % This will create p.local(:).edgeFile & borderFile
    % See graphConstruction subfolder
    job = graphConstruction(p);
    Cluster.waitForJob(job);
    
    % Build segment-to-cube mapping
    job = startCPU( ...
        @Seg.Global.buildSegToCubeMap, ...
        {{p}}, 'segToCubeMap');
    Cluster.waitForJob(job);
    
    % Build segment masks and bounding boxes
    % This produces the 'segMasks.mat' files and is required
    % to produce iso-surfaces for visualizations, etc.
    job = startSegmentMasks(p);
    Cluster.waitForJob(job);
    
    % Calculate edge based features as used in first GP training
    % Calculate segment based features as used in spine head detection
    % see filterbank 
    job = miniFeature(p);
    Cluster.waitForJob(job);
    
    % Make predictions on edge based features using previously trained GP
    job = makePredictions(p,'edges');
    Cluster.waitForJob(job);
    
    % Make predictions for spine heads on segment based features using previously trained spine head classifier
    job = spineHeadDetectionOnCluster(p);
    Cluster.waitForJob(job);
    
    % Run interface classifier using Benedikt's trained classifier and store features
    job = interfaceClassificationOnCluster(p);
    Cluster.waitForJob(job);
    
    % Use interface features to make predictions i.e. generate synapse scores 
    job=makeInterfacePredictionsOnCluster(p);
    Cluster.watiForJob(job);

    %Save the global SVG data
    job = collectSvgDataOnCluster(p);
    Cluster.waitForJob(job);
    
    %Create graph struct 
    job = collectGraphStructOnCluster(p);
    Cluster.waitForJob(job);
    
    

end

