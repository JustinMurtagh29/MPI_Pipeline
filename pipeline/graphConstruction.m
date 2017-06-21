function job = graphConstruction(parameter)
    % collect parameters
    cubeCount = numel(parameter.local);
    inputCell = cell(cubeCount, 1);
    toDelete = false(cubeCount,1);
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
        elseif ~any(mask(:))
            toDelete(curIdx) = 1;
            Util.save(curCube.edgeFile, zeros(0,2,'uint32'));
            Util.save(curCube.borderFile, struct('PixelIdxList',{},'Area',{},'Centroid',{}));
            Util.save(curCube.segmentFile, struct('PixelIdxList',{},'Id',{}));
        else
            inputCell{curIdx} = { ...
                curSegFile, ...
                curCube.edgeFile, ...
                curCube.borderFile, ...
                curCube.segmentFile, ...
                ~mask};
        end
    end
    inputCell(toDelete) = [];
    functionH = @findEdgesAndBordersWrapper;
    job = startCPU(functionH, inputCell, 'graphConstruction', 12); % 30 min should suffice
end

