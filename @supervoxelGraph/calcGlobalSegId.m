
function job = calcGlobalSegId(parameter)

for i=1:size(parameter.local,1)
    for j=1:size(parameter.local,2)
        for k=1:size(parameter.local,3)
            idx = sub2ind(size(parameter.local), i, j, k);
            functionH{idx} = @globalSegId;
            inputCell{idx} = {parameter, i,j,k};
        end
    end
end
job = startCPU(functionH, inputCell, 'GlobalIdCalculation');
end
