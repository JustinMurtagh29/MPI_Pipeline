function splitData = loadSplitData(param, taskGenFile, taskIdFile, nmlDir)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    import connectEM.Chiasma.Util.loadTaskIds;
    import connectEM.Chiasma.Util.loadFlightPaths;
    import connectEM.Chiasma.Flight.buildSplitData;
    
    % task definitions
    taskGenData = load(taskGenFile);
    taskDefs = taskGenData.taskDefs;
    exits = taskGenData.exits;

    % chiasmata
    chiasmata = load(taskGenData.info.param.chiasmataFile);
    axonFile = chiasmata.info.param.axonFile;
    chiasmata = chiasmata.chiasmata;

    % task IDs
    taskIds = loadTaskIds(taskIdFile);
    taskIds = taskIds.id;

    % flight paths
    flights = loadFlightPaths(param, nmlDir);

    % build split file
    splitData = buildSplitData( ...
        chiasmata, taskDefs, exits, taskIds,flights);
    splitData.axonFile = axonFile;
end
