function job = makePredictions(parameter)

for i=1:size(parameter.local,1)
	for j=1:size(parameter.local,2)
		for k=1:size(parameter.local,3)
			if ~exist(parameter.local(i,j,k).saveFolder, 'dir')
				mkdir(parameter.local(i,j,k).saveFolder);
			end
			idx = sub2ind(size(parameter.local), i, j, k);
			functionH{idx} = @edgeProbabilityPrediction;
			inputCell{idx} = {parameter.local(i,j,k).weightFile, parameter.gp.normValues, parameter.gp.initalGroundTruth, parameter.hyperParameter, parameter.local(i,j,k).probFile};
		end
	end
end

job = startCPU(functionH, inputCell, 'graphConstruction');

end

