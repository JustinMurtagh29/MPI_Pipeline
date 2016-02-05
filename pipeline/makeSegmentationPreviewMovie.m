function outputFile = makeSegmentationPreviewMovie(p)
    % Make movie to visualize segmentation given some parameter set p

    % First use same function as bigFwdPass.m to do classification
    functionH = @onlyFwdPass3DonKnossosFolder;
    tempClass.root = [p.tempFolder 'classForMovie/'];
    tempClass.prefix = 'classForMovie';
    inputCell{1} = {p.cnn.first, p.cnn.GPU, p.raw, tempClass, p.bbox};
    if p.cnn.GPU
        job = startGPU(functionH, inputCell, 'classification for movie');
    else
        job = startCPU(functionH, inputCell, 'classification for movie');
    end
    waitForJob(job);
    % Now do segmentation as in miniSegmentation.m
    tempSegFile = [p.tempFolder 'seg.mat'];
	inputCell{1} = {p.class.root p.class.prefix, p.bbox, tempSegFile};
    functionH = parameter.seg.func;
    job = startCPU(functionH, inputCell, 'segmentation for movie');
    waitForJob(job);
    % Now low raw data and segmentation 
    raw = loadRawData(p.raw, p.class, p.bbox, false);
    load(tempSegFile);
    % Make movie
    outputFile = [p.tempFolder datestr(clock, 30)];
    makeSegMovie(seg, raw, outputFile);
end
