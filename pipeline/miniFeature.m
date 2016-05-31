function job = miniFeature(parameter)
% calculates the feature weights and takes 'raw' or 'aff' as input

for i=1:size(parameter.local,1)
	for j=1:size(parameter.local,2)
		for k=1:size(parameter.local,3)
			if ~exist(parameter.local(i,j,k).saveFolder, 'dir')
				mkdir(parameter.local(i,j,k).saveFolder);
            end
            sub = [i,j,k];
			idx = sub2ind(size(parameter.local), i, j, k);
            inputCell{idx} = {parameter, sub, p.norm.func};
		end
	end
end

functionH = parameter.feature.func; 
job = startCPU(functionH, inputCell, 'featureCalculation');

end

