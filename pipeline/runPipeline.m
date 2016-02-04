function runPipeline(p)

    % Runs CNN based forward pass for region defined in p and saves as p.class
    job = bigFwdPass(p);
    waitForJob(job);
    % Runs watershed based segmentation for region defined in p and saves as p.local.segFile
    job = miniSegmentation(p);
    waitForJob(job);

    job = correspondenceFinder(p);

    job = graphConstruction(p);
    waitForJob(job);

end

