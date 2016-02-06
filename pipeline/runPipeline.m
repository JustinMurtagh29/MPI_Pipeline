function runPipeline(p)

    % To be REMOVED!! For 07x2 consistency checks use right normalization values
    p.norm.func = @(x)normalizeStack(x,122,22);
    % Runs CNN based forward pass for region defined in p.bbox,p.raw and saves as p.class
    % Because CNN is translation invariant, saved as KNOSSOS hierachy again
    % Pass bounding box as well as tileSize will be added in bigFwdPass otherwise (legacy stuff)
    job = bigFwdPass(p, p.bbox);
    Cluster.waitForJob(job);
    % Runs watershed based segmentation for region defined in p.bbox,p.raw and saves as p.local.segFile
    % Because watershed segmentation has FOV effects (no translation invariance), processed with large
    % overlap and later joined together
    job = miniSegmentation(p);
    Cluster.waitForJob(job);
    % Find correspondences between tiles processed (segmented) in paralell
    job = correspondenceFinder(p);
    Cluster.waitForJob(job);
    % Make segmentation IDs unique over dataset (were only unique in each
    % tile before), this will be called global IDs from time to time
    job = globalization(p);
    Cluster.waitForJob(job);
    % Construct graph on globalized version of segmentation
    job = graphConstruction(p);
    Cluster.waitForJob(job);
    % Calculate features for edges as used for inital GP    
    job = miniFeature(p);
    Cluster.waitForJob(job);
    
end

