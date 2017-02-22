function job = graphConstruction(parameter)
    % collect parameters
    cubeCount = numel(parameter.local);
    inputCell = cell(cubeCount, 1);

    for curIdx = 1:cubeCount
        curCube = parameter.local(curIdx);
        curSegFile = [curCube.saveFolder, 'segGlobal.mat'];
	if exist(regexprep(parameter.raw.root,'/color/','/mask/'),'dir')  % check for mask data
            mask = readKnossosRoi(regexprep(parameter.raw.root,'/color/','/mask/'),[parameter.raw.prefix,'_mask'],curCube.bboxSmall);
	else
	    mask = 1;
	end
        if all(mask(:))  % ignore mask if everything is one (no padded area)
            inputCell{curIdx} = { ...
                curSegFile, ...
                curCube.edgeFile, ...
                curCube.borderFile, ...
                curCube.segmentFile};
        else
            inputCell{curIdx} = { ...
                curSegFile, ...
                curCube.edgeFile, ...
                curCube.borderFile, ...
                curCube.segmentFile, ...
                ~mask};
        end
    end

    functionH = @findEdgesAndBordersWrapper;
    job = startCPU(functionH, inputCell, 'graphConstruction', 0.4); % 30 min should suffice
end

