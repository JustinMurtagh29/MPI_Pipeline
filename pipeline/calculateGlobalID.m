function job = calcGlobalSegId(parameter)
% calculates number of segments per cube

functionH{1} = @calcNumberSegments;
inputCell{1} = {parameter};

job = startCPU(functionH, inputCell, 'calculateNumberSegments')

for i=1:size(parameter.local,1)
    for j=1:size(parameter.local,2)
        for k=1:size(parameter.local,3)
            idx = sub2ind(size(parameter.local), i, j, k);
            functionH{idx} = @globalSegId;
            inputCell{idx} = {parameter, i,j,k};
        end
    end
end
job = startCPU(functionH, inputCell, 'assignGlobalIDs')

end
