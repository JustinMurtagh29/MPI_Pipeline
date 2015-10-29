function parameterSearchSeg(pT)
    % Sets parameters for segmentation optimization and performs scan for optimal segmentation

    % Output folder for figures and tables
    psValues.outputFolder = [pT.saveFolder 'segmentationPS/'];
    % Radii for Morphological Reconstruction
    psValues.r = [0 1];
    % Parameter for H-minima based segmentation
    psValues.algo(1).fun = @(seg,pars) watershedSeg_v1_cortex(seg, pars(:));
    psValues.algo(1).par = {0:0.05:1 0:50:200};
    % Parameter for threshold based segmentation
    psValues.algo(2).fun = @(seg,pars) watershedSeg_v2_cortex(seg, pars(:));
    psValues.algo(2).par = {-1:0.1:1 0:50:200};
    % from Alex, not sure, check again sometime?
    psValues.algo(3).fun = @(seg,pars) watershedSeg_v3_cortex(seg, pars(:));
    psValues.algo(3).par = {};
    % Same as above, salient watershed idea
    psValues.algo(4).fun = @(seg,pars) salient_watershedSeg_cortex(seg, pars(:));
    psValues.algo(4).par = {};
    % Set parameter for visualization of results
    psValues.makeSegMovies = true;
    psValues.makeErrorMovies = true;
    psValues.makeErrorStacks = true;
    psValues.plotObjSizeHist = true;
    psValues.plotObjChains = true;
    psValues.plotSynRate = true;
    % Calculate all parameter/algo pairs and put into linear structure
    psValues.paramCell = getParamCombinations(psValues.algo);

    if ~exist(psValues.outputFolder, 'dir')
        mkdir(psValues.outputFolder);
    end

    display(['Starting morphological reconstruction:']);
    tic;
    for localIdx=1:length(pT.local)
        for radiusIdx=1:length(psValues.r)
            outputFile = [psValues.outputFolder filesep 'morphRecon-' num2str(localIdx) '-' num2str(radiusIdx)  '.mat'];
            idx = sub2ind([length(pT.local) length(psValues.r)], localIdx, radiusIdx);
            inputCell{idx} = {pT.class, pT.local(localIdx).bboxSmall, psValues.r(radiusIdx), outputFile};
        end
    end

    functionH = @morphologicalReconstruction;
    job = startCPU(functionH, inputCell, 'morphological reconstruction');
    wait(job, 'finished');
    toc;

    clear job functionH inputCell;

    display(['Starting segmentation and evaluation:']);
    tic;
    for localIdx=1:length(pT.local)
        for radiusIdx=1:length(psValues.r)
            for parVar=1:length(psValues.paramCell)
                inputFile = [psValues.outputFolder filesep 'morphRecon-' num2str(localIdx) '-' num2str(radiusIdx)  '.mat'];
                outputFile = [psValues.outputFolder filesep 'seg-' num2str(localIdx) '-' num2str(radiusIdx) '-' num2str(parVar) '.mat'];
                idx = sub2ind([length(pT.local) length(psValues.r) length(psValues.paramCell)], localIdx, radiusIdx, parVar);
                inputCell{idx} = {psValues.paramCell{parVar}, inputFile, outputFile};
            end
        end
    end

    functionH = @performSegmentation;
    job = startCPU(functionH, inputCell, 'segmentation parameter search');
    wait(job, 'finished');
    toc;
    
    display('Parameter search finished. Continuing with visualization of results!');
    tic;
    visualizeParameterSearchSeg(pT, psValues);   
    toc;

end

