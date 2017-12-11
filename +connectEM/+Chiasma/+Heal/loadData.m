function data = loadData(param, taskGenFile, taskIdFile, nmlDir)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    taskGenData = load(taskGenFile);
    taskIds = connectEM.Chiasma.Util.loadTaskIds(taskIdFile);
    flights = connectEM.Flight.loadFromDirs(param, nmlDir);
    
    data = struct;
    data.axonFile = taskGenData.info.param.axonFile;
    data.endings  = taskGenData.endings;
    data.taskDef  = taskGenData.taskDefs;
    data.taskIds  = taskIds.id;
    data.flights  = flights;
end