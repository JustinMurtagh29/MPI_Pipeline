function job = makePredictions(parameter,mode)

%copy GP state variables to results/state
if ~exist(parameter.gp.hyperParameter,'file') || ~exist(parameter.gp.initalGroundTruth,'file')
me = mfilename;
mydir = which(me);
mydir = mydir(1:end-2-numel(me)); 
copyfile([mydir '../state/'],parameter.gp.stateFolder);
end

%Generate normValues.mat based on quantiles
calculateNormValues(parameter);

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
                inputCell{idx} = {parameter.local(i,j,k).weightFile, parameter.gp.normValues, parameter.gp.initialGroundTruth, parameter.gp.hyperParameter, parameter.local(i,j,k).probFile};
            elseif strcmp(mode,'glia')  
                functionH{idx} = @edgeProbabilityPrediction;
                inputCell{idx} = {parameter.local(i,j,k).segmentWeightFile, parameter.glia.normValues, parameter.glia.initialGroundTruth, parameter.glia.hyperParameter, parameter.local(i,j,k).gliaProbFile};
            end			
		end
	end
end

job = startCPU(functionH, inputCell, 'makePredictions');

end

