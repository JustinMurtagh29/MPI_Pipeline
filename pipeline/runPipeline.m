function runPipeline(p)

    % To be REMOVED!! For 07x2 consistency checks use right normalization values
    p.norm.func = @(x)normalizeStack(x,122,22);
    % Runs CNN based forward pass for region defined in p.bbox,p.raw and saves as p.class
    % Because CNN is translation invariant, saved as KNOSSOS hierachy again
    % Pass bounding box as well as tileSize will be added in bigFwdPass otherwise (legacy stuff)
    % This uses CNN subfolder in code repository
    job = bigFwdPass(p, p.bbox);
    Cluster.waitForJob(job);
    % Runs watershed based segmentation for region defined in p.bbox on p.raw and saves as p.local.segFile
    % Because watershed segmentation has FOV effects (no translation invariance), processed with large
    % overlap and later joined together (see correspondences)
    % Uses segmentation subfolder in code repository
    job = miniSegmentation(p);
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
    jobs = globalization(p);
    Cluster.waitForJob(jobs);
    % Construct graph on globalized version of segmentation
    % This will create p.local(:).edgeFile & borderFile
    % See graphConstruction subfolder
    job = graphConstruction(p);
    Cluster.waitForJob(job);
    % Calculate features for edges as used for inital GP
    % Calculate edge based features as used in first GP training
    % see filterbank 
    job = miniFeature(p);
    Cluster.waitForJob(job);
    
end

