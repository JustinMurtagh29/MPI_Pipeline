function runPipeline(p)

    % Runs CNN based forward pass for region defined in p.bbox,p.raw and saves as p.class
    % Because CNN is translation invariant, saved as KNOSSOS hierachy again
    job = bigFwdPass(p);
    waitForJob(job);
    % Runs watershed based segmentation for region defined in p.bbox,p.raw and saves as p.local.segFile
    % Because watershed segmentation has FOV effects (no translation invariance), processed with large
    % overlap and later joined together
    job = miniSegmentation(p);
    waitForJob(job);
    % Find correspondences between tiles processed (segmented) in paralell
    job = correspondenceFinder(p);
    waitForJob(job);
    % Make segmentation IDs unique over dataset (were only unique in each
    % tile before), this will be called global IDs from time to time
    job = globalization(p);
    waitForJob(job);
    % Construct graph on globalized version of segmentation
    job = graphConstruction(p);
    waitForJob(job);
    % Calculate features for edges as used for inital GP    
    job = miniFeature(p);
    waitForJob(job);
    
end

