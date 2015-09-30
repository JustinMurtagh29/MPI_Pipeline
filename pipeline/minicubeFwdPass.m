function job = minicubeFwdPass( parameter )
    % classification of subergions within the data set (for segmentation optimization and GP training) 

    if isfield(parameter.cnn, 'benedikt') && parameter.cnn.benedikt;
        load(parameter.cnn.saveFile);
        for tr=1:length(parameter.local)
            bbox = parameter.local(tr).bboxBig;
            options.gpuDev = parameter.cnn.GPU;
            options.val_fwd_alg = 'memory1';
            options.target_size = [128 128 128];
            inputCell{tr} = {cnet, bbox, options, parameter.raw, parameter.local(tr).class};
        end
        functionH = @Codat.CNN.Misc.predictROI_forCluster; 
    else
        for tr=1:length(parameter.local)
            bbox = parameter.local(tr).bboxBig;
            inputCell{tr} = {parameter.cnn.first, parameter.cnn.GPU, parameter.raw, parameter.local(tr).class, bbox};
        end
        functionH = @onlyFwdPass3DonKnossosFolder;
    end

    if parameter.cnn.GPU
        job = startGPU(functionH, inputCell, 'classification');
    else
        job = startCPU(functionH, inputCell, 'classification');
    end

end

