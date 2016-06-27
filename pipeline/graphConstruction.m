function job = graphConstruction(parameter)
    % collect parameters
    cubeCount = numel(parameter.local);
    inputCell = cell(cubeCount, 1);

    for curIdx = 1:cubeCount
        curCube = parameter.local(curIdx);
        curSegFile = [curCube.saveFolder, 'segGlobal.mat'];

        inputCell{curIdx} = { ...
            curSegFile, ...
            curCube.edgeFile, ...
            curCube.borderFile, ...
            curCube.segmentFile};
    end

    functionH = @findEdgesAndBordersFast;
    job = startCPU(functionH, inputCell, 'graphConstruction');
end

