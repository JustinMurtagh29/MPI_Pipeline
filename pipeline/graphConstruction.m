function job = graphConstruction(parameter)

for i=1:size(parameter.local,1)
	for j=1:size(parameter.local,2)
		for k=1:size(parameter.local,3)
			if ~exist(parameter.local(i,j,k).saveFolder, 'dir')
				mkdir(parameter.local(i,j,k).saveFolder);
			end
			idx = sub2ind(size(parameter.local), i, j, k);
			functionH{idx} = @findEdgesAndBordersFast; 
			inputCell{idx} = {parameter.local(i,j,k).segFile, parameter.local(i,j,k).edgeFile, parameter.local(i,j,k).borderFile, parameter.local(i,j,k).segmentFile, parameter.tileBorder};
		end
	end
end

job = startCPU(functionH, inputCell, 'graphConstruction');

end

