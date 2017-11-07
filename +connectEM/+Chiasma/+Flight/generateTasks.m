function taskDef = generateTasks(param, chiasmata, exits, outputFile)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    %
    % Based on code by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Kevin M. Boergens <kevin.boergens@brain.mpg.de>
    
    %% configuration
    boxMarginNm = 200;
    voxelSize = param.raw.voxelSize;
        
    %% main
    exitCount = size(exits, 1);
    
    taskDef = table;
    taskDef.position = zeros(exitCount, 3);
    taskDef.rotation = zeros(exitCount, 3);
    taskDef.bbox = zeros(exitCount, 6);
    
    for curIdx = 1:exitCount
        curExit = exits(curIdx, :);
        curChiasma = chiasmata{curExit.aggloId};
        
        curDirections = curChiasma.direction{curExit.chiasmaId};
        curPositions = curChiasma.position{curExit.chiasmaId};
        curPosition = curPositions(curExit.exitId, :);
        
        % determine Euler angles
        curEulerAngles = calculateEulerAngles( ...
            -curDirections(curExit.exitId, :), voxelSize);
        
        % find what bounding box to use
        curMaxDistNm = pdist2( ...
            curPosition .* voxelSize, ...
            curPositions .* voxelSize, ...
            'squaredeuclidean');
        curMaxDistNm = sqrt(max(curMaxDistNm));
        
        curBox = ceil((curMaxDistNm + boxMarginNm) ./ voxelSize);
        curBox = curPosition(:) + ([-1, +1] .* curBox(:));
        
        % build webKNOSSOS-style bounding box
        curBoxWk = cat(2, curBox(:, 1)', 1 + diff(curBox, 1, 2)');
        
        taskDef.position(curIdx, :) = curPosition - 1;
        taskDef.rotation(curIdx, :) = curEulerAngles;
        taskDef.bbox(curIdx, :)     = curBoxWk;
    end
    
    %% build output
    taskParam = struct;
    taskParam.dataSet = '2012-09-28_ex145_07x2_ROI2017';
    taskParam.taskTypeId = '56d6a7c6140000d81030701e';
    taskParam.expDomain = 'queriesMHchiasma';
    taskParam.expMinVal = 1;
    taskParam.instances = 1;
    taskParam.team = 'Connectomics department';
    taskParam.project = 'queriesMHchiasma';
    
    % write webKNOSSOS task definitions
    taskDef = ...
        connectEM.exportTaskDefinitions(taskParam, taskDef, outputFile);
end

function eulerAngles = calculateEulerAngles(di, voxelSize)
    % Make sure direction is normalized
    di = di .* voxelSize;
    di = di ./ norm(di);
    
    % Calculate euler angles, from wk-paper -> RESCOPaaS
    eulerAngles = connectEM.diffToEulerAngle(di);
    eulerAngles = round(reshape(eulerAngles, 1, []));
end