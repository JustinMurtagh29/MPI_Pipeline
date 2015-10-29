function job = minicubeFwdPass( parameter )
    % classification of subergions within the data set (for segmentation optimization and GP training) 

    if isfield(parameter.cnn, 'benedikt') && parameter.cnn.benedikt;
        load(parameter.cnn.saveFile);
        % Settings for Benedikt's CNN
        options.gpuDev = parameter.cnn.GPU;
        options.val_fwd_alg = 'memory1';
        options.target_size = [128 128 128];
        % Initalize counter variable (each KNOSSOS cube will be one task)
        nrTasks = 0;
        for tr=1:length(parameter.local)
            bbox = parameter.local(tr).bboxBig;
            % Bbox inidcating lower,left,front corner of Knossos cubes
            bbox = floor((bbox-1)./128).*128+1;
            for x=bbox(1,1):128:bbox(1,2)
                for y=bbox(2,1):128:bbox(2,2)
                    for z=bbox(3,1):128:bbox(3,2)
                        % Bbox of this Knossos cube
                        thisBbox = [x x+127; y y+127; z z+127];
                        nrTasks = nrTasks + 1;
                        inputCell{nrTasks} = {cnet, thisBbox, options, parameter.raw, parameter.local(tr).class};
                    end
                end
            end
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

