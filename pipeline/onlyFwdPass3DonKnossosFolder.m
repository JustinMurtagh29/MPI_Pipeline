function onlyFwdPass3DonKnossosFolder(cnetLocation, gpuSwitch, input, result, bbox)

% Choose the right GPU if on actvtClass indicated GPU usage
w = getCurrentWorker();
w = w.Name;
gpuDevice(str2double(w(14:15)));

% Load the CNN used for classification
load(cnetLocation, 'cnet');
cnet = cnet.loadLastCNN;
if gpuSwitch
    cnet.run.actvtClass = @gpuArray;
else
    cnet.run.actvtClass = @single;
end

% Load data with right border for cnet
bboxWithBorder(:,1) = bbox(:,1) - ceil(cnet.randOfConvn'/2);
bboxWithBorder(:,2) = bbox(:,2) + ceil(cnet.randOfConvn'/2);
raw = readKnossosRoi(input.root, input.prefix, bboxWithBorder);

% Normalize data
if cnet.normalize
	raw = normalizeStack(single(raw));
else
	raw = single(raw);
end

% Memory efficent fwd pass
classification = onlyFwdPass3D(cnet, raw);
% Save result to KNOSSOS folder
writeKnossosRoi(result.root, result.prefix, bbox(:,1)', classification, 'single');

end
