function taskDef = generateTasks(param, axons, outputFile)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % build task definitions
    taskDef = connectEM.Chiasma.Heal.buildTaskDefs(param, axons);
    
    %% set bounding box
    boxOff = param.bbox(:, 1)';
    boxSize = 1 + diff(param.bbox, 1, 2)';
    bbox = horzcat(boxOff, boxSize);
    
    % restrict flights to segmentation
    taskDef.bbox = repmat(bbox, size(taskDef, 1), 1);

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