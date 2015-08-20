function job = minicubeFwdPass( parameter )
    % classification of subergions within the data set (for segmentation optimization and GP training) 

    if isfield(cnn, 'benedikt') && cnn.benedikt
        import benedikt;
        load parameter.cnn.saveFile;
        for tr=1:length(parameter.local)
            bbox = parameter.local(tr).bboxBig;
            bbox(:,1) - mod(bbox(:,1)-1,128);
            bbox(:,2) = bbox(:,2) + (128 - mod(bbox(:,2),128));

            inputCell{tr} = {parameter.n.savefirst, parameter.cnn.GPU, parameter.raw, parameter.local(tr).class, bbox};
        end
        functionH = @benedikt.Codat.CNN.Misc.predictROI; 
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

