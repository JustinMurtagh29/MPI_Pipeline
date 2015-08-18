function job = minicubeFwdPass( parameter )
% classification of subergions within the data set (for segmentation optimization and GP training) 

for tr=1:length(parameter.local)
	bbox = parameter.local(tr).bboxBig;
	bbox(:,1) - mod(bbox(:,1)-1,128);
	bbox(:,2) = bbox(:,2) + (128 - mod(bbox(:,2),128));
    inputCell{tr} = {parameter.cnn.first, parameter.cnn.GPU, parameter.raw, parameter.local(tr).class, bbox};
end

functionH = @onlyFwdPass3DonKnossosFolder;

if parameter.cnn.GPU
	job = startGPU(functionH, inputCell, 'classification');
else
	job = startCPU(functionH, inputCell, 'classification');
end

end

