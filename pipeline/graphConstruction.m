function job = graphConstruction(parameter)

for i=1:size(parameter.local,1)
	for j=1:size(parameter.local,2)
		for k=1:size(parameter.local,3)
			idx = sub2ind(size(parameter.local), i, j, k);
			inputCell{idx} = {[parameter.local(i,j,k).saveFolder 'segGlobal.mat'], ...
                parameter.local(i,j,k).edgeFile, parameter.local(i,j,k).borderFile, ...
                parameter.local(i,j,k).segmentFile};
		end
	end
end

functionH = @findEdgesAndBordersFast; 
job = startCPU(functionH, inputCell, 'graphConstruction');

end

