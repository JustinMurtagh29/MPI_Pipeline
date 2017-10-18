function exportTaskDefinitions(param, taskDef, outFile)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    fullTaskDef = taskDef;
    fullTaskDef.dataSet(:) = {param.dataset};
    fullTaskDef.taskTypeId(:) = {param.taskTypeId};
    fullTaskDef.expDomain(:) = {param.expDomain};
    fullTaskDef.expMinVal(:) = param.expMinVal;
    fullTaskDef.instances(:) = param.instances;
    fullTaskDef.team(:) = {param.team};
    fullTaskDef.project(:) = {param.project};
    
    % reorder-columns
    columnOrder = { ...
        'dataSet', 'taskTypeId', 'expDomain', 'expMinVal', ...
        'position', 'rotation', 'instances', 'team', 'bbox', 'project'};
    fullTaskDef = fullTaskDef(:, columnOrder);
    
    % write CSV file
    writetable(fullTaskDef, outFile, 'WriteVariableNames', false);
end