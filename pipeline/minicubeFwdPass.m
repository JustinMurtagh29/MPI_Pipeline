function job = minicubeFwdPass( parameter )
    % classification of subergions within the data set (for segmentation optimization and GP training) 

    if isfield(parameter.cnn, 'benedikt') && parameter.cnn.benedikt;
        load(parameter.cnn.saveFile);
        for tr=1:length(parameter.local)
            bbox = parameter.local(tr).bboxBig;
            % Cant remember why this should be needed. Drop?
            bbox(:,1) - mod(bbox(:,1)-1,128);
            bbox(:,2) = bbox(:,2) + (128 - mod(bbox(:,2),128));
            options.gpuDev = parameter.cnn.GPU;
            options.target_size = [128 128 128];
            inputCell{tr} = {cnet, bbox, options, parameter.raw, parameter.local(tr).class};
        end
        functionH = @Codat.CNN.Misc.predictROI_forCluster; 
    else
        for tr=1:length(parameter.local)
            bbox = parameter.local(tr).bboxBig;
            bbox(:,1) - mod(bbox(:,1)-1,128);
            bbox(:,2) = bbox(:,2) + (128 - mod(bbox(:,2),128));
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

