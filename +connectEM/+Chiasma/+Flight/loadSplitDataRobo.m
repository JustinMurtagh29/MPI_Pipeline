function splitData = loadSplitDataRobo(param, taskGenFile, roboRunIds)
    % Written by
    %   Martin Schmidt <martin.schmidt@brain.mpg.de>
    import connectEM.RoboFlight.loadFromRoboTDefs;
    import connectEM.Chiasma.Flight.buildSplitData;
    
    % task definitions
    taskGenData = load(taskGenFile);
    taskDefs = taskGenData.taskDefs;
    exits = taskGenData.exits;
    roboem = taskGenData.roboem;

    % chiasmata
    chiasmataFile = taskGenData.chiasmataFile;
    chiasmata = load(chiasmataFile);
    axonFile = chiasmata.info.param.axonFile;
    chiasmata = chiasmata.chiasmata;

    % flight paths
    flights = loadFromRoboTDefs(roboem.tdefOutputDir, roboem.tdefInfo.batch_params.tdefFileNames, roboRunIds, param);

    % build split file
    splitData = buildSplitData( ...
        chiasmata, taskDefs, exits, flights.globTaskIds, flights);
    splitData.axonFile = axonFile;
end

