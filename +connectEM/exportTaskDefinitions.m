function fullTaskDef = exportTaskDefinitions(param, taskDef, outFile)
    % fullTaskDef = exportTaskDefinitions(param, taskDef, outFile)
    %   Completes the flight task definitions in `taskDef` with the default
    %   values in `param`, and writes the result to `outFile`.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    fullTaskDef = taskDef;
    fullTaskDef.dataSet(:) = {param.dataSet};
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
    
    if exist('outFile', 'var') && ~isempty(outFile)
        % write results to CSV file, if file path is specified
        writetable(fullTaskDef, outFile, 'WriteVariableNames', false);
    end
end