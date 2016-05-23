function runPipeline(p)

    % Runs CNN based forward pass for region defined in p.bbox,p.raw and saves as p.class	
    % Because CNN is translation invariant, saved as KNOSSOS hierachy again 
    % Pass bounding box as well as tileSize will be added in bigFwdPass otherwise (legacy stuff)
    % This uses CNN subfolder in code repository
    if p.retina
    job= bigFwdPassRetina(p,p.bbox);
    else
    job = bigFwdPass(p, p.bbox);
    end
    Cluster.waitForJob(job);
    
    % Runs watershed based segmentation for region defined in p.bbox on p.raw and saves as p.local.segFile
    % Because watershed segmentation has FOV effects (no translation invariance), processed with large
    % overlap and later joined together (see correspondences)
    % Uses segmentation subfolder in code repository
    if p.retina
    job = miniSegmentationRetina(p);
    else   
    job = miniSegmentation(p);
    end
    Cluster.waitForJob(job);
    
    % Find correspondences between tiles processed (segmented) in paralell
    % Uses correspondences subfolder in code repository
    job = correspondenceFinder(p);
    Cluster.waitForJob(job);
    
    % Transfer segmentation from tempSegFile to segFile and drop overlaps 
    % Also in correspondences subfolder
    job = removeOverlaps(p);
    Cluster.waitForJob(job);
    
    % Make segmentation IDs unique over dataset (were only unique in each
    % tile before), this will be called global IDs from time to time
    % will work on correspondences and segmentation
    % see globalization subfolder 
    % globalization is run in two parts: globalizeSegmentation and globalizeCorrespondences
    job = globalizeSegmentation(p);
    Cluster.waitForJob(job);
    
    job = globalizeCorrespondences(p);
    Cluster.waitForJob(job);
    
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
    
    % Calculate features for edges as used for inital GP
    % Calculate edge based features as used in first GP training
    % see filterbank 
    job = miniFeature(p);
    Cluster.waitForJob(job);
    
    %Make predictions on edge based features using previously trained GP
    job = makePredictions(p,'edges');
    Cluster.waitForJob(job);
    
    %Run interface classifier using Benedikt's trained classifier
    job = interfaceClassificationOnCluster(p);
    Cluster.waitForJob(job);

    %Save the global SVG data
    job = collectSvgDataOnCluster(p);
    Cluster.waitForJob(job);
    
    %Create graph struct 
    job = collectGraphStructOnCluster(p);
    Cluster.waitForJob(job);
    
    

end

