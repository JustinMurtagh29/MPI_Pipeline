function splitData = loadSplitData(param, taskGenFile, taskIdFile, nmlDir)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    import connectEM.Chiasma.Util.loadTaskIds;
    import connectEM.Chiasma.Flight.buildSplitData;
    import connectEM.Flight.loadFromDirs;
    
    % task definitions
    taskGenData = load(taskGenFile);
    taskDefs = taskGenData.taskDefs;
    exits = taskGenData.exits;

    % chiasmata
    chiasmataFile = taskGenData.chiasmataFile;
    chiasmata = load(chiasmataFile);
    axonFile = chiasmata.info.param.axonFile;
    chiasmata = chiasmata.chiasmata;

    % task IDs
    taskIds = loadTaskIds(taskIdFile);
    taskIds = taskIds.id;

    % flight paths
    flights = loadFromDirs(param, nmlDir);

    % build split file
    splitData = buildSplitData( ...
        chiasmata, taskDefs, exits, taskIds, flights);
    splitData.axonFile = axonFile;
end
