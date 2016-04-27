function job = makePredictions(parameter,mode)


% Add visualization of some statistics

for i=1:size(parameter.local,1)
	for j=1:size(parameter.local,2)
		for k=1:size(parameter.local,3)
			if ~exist(parameter.local(i,j,k).saveFolder, 'dir')
				mkdir(parameter.local(i,j,k).saveFolder);
			end
			idx = sub2ind(size(parameter.local), i, j, k);
            if strcmp(mode,'edges')
                functionH{idx} = @edgeProbabilityPrediction;
                inputCell{idx} = {parameter.local(i,j,k).weightFile, parameter.gp.normValues, parameter.gp.initalGroundTruth, parameter.gp.hyperParameter, parameter.local(i,j,k).probFile};
            elseif strcmp(mode,'glia')  
                functionH{idx} = @edgeProbabilityPrediction;
                inputCell{idx} = {parameter.local(i,j,k).segmentWeightFile, parameter.glia.normValues, parameter.glia.initalGroundTruth, parameter.glia.hyperParameter, parameter.local(i,j,k).gliaProbFile};
            end			
		end
	end
end

job = startCPU(functionH, inputCell, 'graphConstruction');

end

