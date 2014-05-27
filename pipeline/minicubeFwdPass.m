function minicubeFwdPass( parameter )
% classification of subergions within the data set (for segmentation optimization and GP training) 

load(parameter.cnn.first, 'cnet');
cnet = cnet.loadLastCNN;
if parameter.cnn.GPU
	cnet.run.actvtClass = @gsingle;
else
	cnet.run.actvtClass = @single;
end

idx = 1;
for tr=1:length(parameter.local)
	bbox = parameter.local(tr).bboxBig;
	bbox(:,1) - mod(bbox(:,1)-1,128);
	bbox(:,2) = bbox(:,2) + (128 - mod(bbox(:,2),128));
	cubeSize = [128 128 128]; % This is not of any importance due to CNN translation invariance
	X = bbox(1,1):cubeSize(1):bbox(1,2)+1;
	Y = bbox(2,1):cubeSize(2):bbox(2,2)+1;
	Z = bbox(3,1):cubeSize(3):bbox(3,2)+1;
	for i=1:length(X)-1
	    for j=1:length(Y)-1
	        for k=1:length(Z)-1
	            functionH{idx} = @onlyFwdPass3DonKnossosFolder;
	            inputCell{idx} = {cnet, parameter.raw, parameter.local(tr).class, [X(i) X(i+1)-1; Y(j) Y(j+1)-1; Z(k) Z(k+1)-1]};
		    idx = idx + 1;
	        end
	    end
	end
end

if parameter.cnn.GPU
	startGPU(functionH, inputCell, 'classification');
else
	startCPU(functionH, inputCell, 'classification');
end

end

